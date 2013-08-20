/*
 * loadGraph.c
 *
 * Copyright (c) 2008-2012 BGI-Shenzhen <soap at genomics dot org dot cn>.
 *
 * This file is part of SOAPdenovo.
 *
 * SOAPdenovo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SOAPdenovo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SOAPdenovo.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "stdinc.h"
#include "newhash.h"
#include "kmerhash.h"
#include "extfunc.h"
#include "extvab.h"

#define preARCBLOCKSIZE 100000

static unsigned int loadArcs ( char * graphfile );
static void loadContig ( char * graphfile );
static int maskRepeatByArc ( unsigned avg_weight );

/*
void loadUpdatedVertex (char *graphfile)
{
    char name[256], line[256];
    FILE *fp;
    Kmer word, bal_word;
    int num_kmer, i;
    char ch;

    sprintf (name, "%s.updated.vertex", graphfile);
    fp = ckopen (name, "r");

    while (fgets (line, sizeof (line), fp) != NULL)
    {
        if (line[0] == 'V')
        {
            sscanf (line + 6, "%d %c %d", &num_kmer, &ch, &overlaplen);
            printf ("there're %d kmers in vertex file\n", num_kmer);
            break;
        }
    }

    vt_array = (VERTEX *) ckalloc ((2 * num_kmer) * sizeof (VERTEX));

    for (i = 0; i < num_kmer; i++)
    {
        fscanf (fp, "%llx ", &word);
        vt_array[i].kmer = word;
    }

    fclose (fp);

    for (i = 0; i < num_kmer; i++)
    {
        bal_word = reverseComplement (vt_array[i].kmer, overlaplen);
        vt_array[i + num_kmer].kmer = bal_word;
    }

    num_vt = num_kmer;
}*/

int uniqueLenSearch ( unsigned int * len_array, unsigned int * flag_array, int num, unsigned int target )
{
	int mid, low, high;
	low = 1;
	high = num;

	while ( low <= high )
	{
		mid = ( low + high ) / 2;

		if ( len_array[mid] == target )
		{
			break;
		}
		else if ( target > len_array[mid] )
		{
			low = mid + 1;
		}
		else
		{
			high = mid - 1;
		}
	}

	if ( low > high )
	{
		return -1;
	}

	//locate the first same length unflaged
	return flag_array[mid]++;
}

int lengthSearch ( unsigned int * len_array, unsigned int * flag_array, int num, unsigned int target )
{
	int mid, low, high, i;
	low = 1;
	high = num;

	while ( low <= high )
	{
		mid = ( low + high ) / 2;

		if ( len_array[mid] == target )
		{
			break;
		}
		else if ( target > len_array[mid] )
		{
			low = mid + 1;
		}
		else
		{
			high = mid - 1;
		}
	}

	if ( low > high )
	{
		return -1;
	}

	//locate the first same length unflaged
	if ( !flag_array[mid] )
	{
		for ( i = mid - 1; i > 0; i-- )
		{
			if ( len_array[i] != len_array[mid] || flag_array[i] )
			{
				break;
			}
		}

		flag_array[i + 1] = 1;
		return i + 1;
	}
	else
	{
		for ( i = mid + 1; i <= num; i++ )
		{
			if ( !flag_array[i] )
			{
				break;
			}
		}

		flag_array[i] = 1;
		return i;
	}
}

void quick_sort_int ( unsigned int * length_array, int low, int high )
{
	int i, j;
	unsigned int pivot;

	if ( low < high )
	{
		pivot = length_array[low];
		i = low;
		j = high;

		while ( i < j )
		{
			while ( i < j && length_array[j] >= pivot )
			{
				j--;
			}

			if ( i < j )
			{
				length_array[i++] = length_array[j];
			}

			while ( i < j && length_array[i] <= pivot )
			{
				i++;
			}

			if ( i < j )
			{
				length_array[j--] = length_array[i];
			}
		}

		length_array[i] = pivot;
		quick_sort_int ( length_array, low, i - 1 );
		quick_sort_int ( length_array, i + 1, high );
	}
}

static int maskRepeatByArc ( unsigned avg_weight )
{
	unsigned int i, bal_i;
	int counter = 0;
	int arc_num;
	unsigned int arc_weight1, arc_weight2;
	preARC * arc;

	for ( i = 1; i <= num_ctg; ++i )
	{
		if ( contig_array[i].mask == 1 )
		{
			if ( isSmallerThanTwin ( i ) )
			{
				++i;
			}

			continue;
		}

		bal_i = getTwinCtg ( i );
		arc = contig_array[bal_i].arcs;
		arc_weight1 = maxArcWeight ( arc );
		arc = contig_array[i].arcs;
		arc_weight2 = maxArcWeight ( arc );

		if ( arc_weight1 + arc_weight2 >= 4 * avg_weight )
		{
			contig_array[i].mask = 1;
			contig_array[bal_i].mask = 1;

			if ( i == bal_i ) { counter += 1; }
			else { counter += 2; }
		}

		if ( isSmallerThanTwin ( i ) )
		{
			++i;
		}
	}

	return counter;
}

/*************************************************
 Function:
    loadUpdatedEdges
 Description:
    Loads contig information and masks some contigs according to setting.
 Input:
    1. graphfile:       prefix of graph file
 Output:
    None.
 Return:
    None.
 *************************************************/
void loadUpdatedEdges ( char * graphfile )
{
	char c, name[256], line[1024];
	int bal_ed, cvg;
	FILE * fp, *out_fp;
	Kmer from_kmer, to_kmer;
	unsigned int num_ctgge, length, index = 0, num_kmer;
	unsigned int i = 0, j;
	int newIndex;
	unsigned int * length_array, *flag_array, diff_len;
	char * outfile = graphfile;
	long long cvgSum = 0;
	long long counter = 0;
	unsigned int avg_arc_wt;
	int ctg_short_cutoff;
	float high_cvg_cutoff1, high_cvg_cutoff2, low_cvg_cutoff;
	int cut_len;
	//get overlaplen from *.preGraphBasic
	sprintf ( name, "%s.preGraphBasic", graphfile );
	fp = ckopen ( name, "r" );

	while ( fgets ( line, sizeof ( line ), fp ) != NULL )
	{
		if ( line[0] == 'V' )
		{
			sscanf ( line + 6, "%d %c %d", &num_kmer, &c, &overlaplen );
			fprintf ( stderr, "Kmer size: %d\n", overlaplen );
			break;
		}
	}

	cut_len = COMPATIBLE_MODE == 0 ? overlaplen : 0;

	if ( ctg_short == 0 )
	{
		ctg_short = overlaplen + 2;
	}

	ctg_short_cutoff = 2 * overlaplen + 2 < 100 ? 100 : 0;
	fclose ( fp );
	sprintf ( name, "%s.updated.edge", graphfile );
	fp = ckopen ( name, "r" );
	sprintf ( name, "%s.newContigIndex", outfile );
	out_fp = ckopen ( name, "w" );

	while ( fgets ( line, sizeof ( line ), fp ) != NULL )
	{
		if ( line[0] == 'E' )
		{
			sscanf ( line + 5, "%d", &num_ctgge );
			fprintf ( stderr, "There are %d edge(s) in edge file.\n", num_ctgge );
			break;
		}
	}

	index_array = ( unsigned int * ) ckalloc ( ( num_ctgge + 1 ) * sizeof ( unsigned int ) );
	length_array = ( unsigned int * ) ckalloc ( ( num_ctgge + 1 ) * sizeof ( unsigned int ) );
	flag_array = ( unsigned int * ) ckalloc ( ( num_ctgge + 1 ) * sizeof ( unsigned int ) );

	while ( fgets ( line, sizeof ( line ), fp ) != NULL )
	{
		if ( line[0] == '>' )
		{
			sscanf ( line + 7, "%d", &length );
			index_array[++index] = length;
			length_array[++i] = length;
		}
	}

	num_ctg = index;
	orig2new = 1;
	qsort ( & ( length_array[1] ), num_ctg, sizeof ( length_array[0] ), cmp_int );
	//extract unique length
	diff_len = 0;

	for ( i = 1; i <= num_ctg; i++ )
	{
		for ( j = i + 1; j <= num_ctg; j++ )
			if ( length_array[j] != length_array[i] )
			{
				break;
			}

		length_array[++diff_len] = length_array[i];
		flag_array[diff_len] = i;
		i = j - 1;
	}

	contig_array = ( CONTIG * ) ckalloc ( ( num_ctg + 1 ) * sizeof ( CONTIG ) );
	//load edges
	index = 0;
	rewind ( fp );

	while ( fgets ( line, sizeof ( line ), fp ) != NULL )
	{
		if ( line[0] == '>' )
		{
			sscanf ( line, ">length %u,%d,%d", &length, &bal_ed, &cvg );
			newIndex = uniqueLenSearch ( length_array, flag_array, diff_len, length );
			index_array[++index] = newIndex;

			if ( length != 0 ) { contig_array[newIndex].length = length - cut_len; }
			else  { contig_array[newIndex].length = 0; }

			contig_array[newIndex].bal_edge = bal_ed + 1;
			contig_array[newIndex].downwardConnect = NULL;
			contig_array[newIndex].mask = 0;
			contig_array[newIndex].flag = 0;
			contig_array[newIndex].arcs = NULL;
			contig_array[newIndex].seq = NULL;
			contig_array[newIndex].multi = 0;
			contig_array[newIndex].inSubGraph = 0;
			contig_array[newIndex].bubbleInScaff = 0;
			contig_array[newIndex].cvg = cvg / 10;

			if ( cvg && length > 100 )
			{
				counter += length - cut_len;
				cvgSum += cvg * ( length - cut_len );
			}

			fprintf ( out_fp, "%d %d %d\n", index, newIndex, contig_array[newIndex].bal_edge );
		}
	}

	if ( counter )
	{
		cvgAvg = cvgSum / counter / 10 > 2 ? cvgSum / counter / 10 : 3;
	}

	//mark repeats
	int bal_i;

	if ( maskRep )
	{
		high_cvg_cutoff1 = cvg_high * cvgAvg;
		high_cvg_cutoff2 = cvg_high * cvgAvg * 0.8;
		low_cvg_cutoff = cvg_low * cvgAvg;
		counter = 0;
		fprintf ( stderr, "Mask contigs with coverage lower than %.1f or higher than %.1f, and strict length %d.\n", low_cvg_cutoff, high_cvg_cutoff1, ctg_short_cutoff );

		for ( i = 1; i <= num_ctg; i++ )
		{
			bal_i = getTwinCtg ( i );

			if ( ( contig_array[i].cvg + contig_array[bal_i].cvg ) > 2 * high_cvg_cutoff1 )
			{
				contig_array[i].mask = 1;
				contig_array[bal_i].mask = 1;

				if ( i == bal_i ) { counter += 1; }
				else { counter += 2; }
			}
			else if ( contig_array[i].length < ctg_short_cutoff && ( contig_array[i].cvg > high_cvg_cutoff2 || contig_array[bal_i].cvg > high_cvg_cutoff2 || ( contig_array[i].cvg < low_cvg_cutoff && contig_array[bal_i].cvg < low_cvg_cutoff ) ) )
			{
				contig_array[i].mask = 1;
				contig_array[bal_i].mask = 1;

				if ( i == bal_i ) { counter += 1; }
				else { counter += 2; }
			}
			else if ( cvgAvg < 50 && ( contig_array[i].cvg >= 63 || contig_array[bal_i].cvg >= 63 ) )
			{
				contig_array[i].mask = 1;
				contig_array[bal_i].mask = 1;

				if ( i == bal_i ) { counter += 1; }
				else { counter += 2; }
			}

			if ( isSmallerThanTwin ( i ) )
			{
				i++;
			}
		}

		fprintf ( stderr, "Average contig coverage is %d, %lld contig(s) masked.\n", cvgAvg, counter );
	}

	counter = 0;

	for ( i = 1; i <= num_ctg; i++ )
	{
		if ( contig_array[i].mask )
		{
			continue;
		}

		bal_i = getTwinCtg ( i );

		if ( contig_array[i].length < ctg_short )
		{
			contig_array[i].mask = 1;
			contig_array[bal_i].mask = 1;

			if ( i == bal_i ) { counter += 1; }
			else { counter += 2; }
		}

		if ( isSmallerThanTwin ( i ) )
		{
			i++;
		}
	}

	fprintf ( stderr, "Mask contigs shorter than %d, %lld contig(s) masked.\n", ctg_short, counter );
	avg_arc_wt = loadArcs ( graphfile );
	counter = 0;
	//counter = maskRepeatByArc(avg_arc_wt);
	//printf ("Mask contigs with multi arcs, %d contig masked\n", counter);
	//tipsCount();
	loadContig ( graphfile );
	fprintf ( stderr, "Done loading updated edges.\n" );
	free ( ( void * ) length_array );
	free ( ( void * ) flag_array );
	fclose ( fp );
	fclose ( out_fp );
}

static void add1Arc ( unsigned int from_ed, unsigned int to_ed, unsigned int weight )
{
	preARC * parc;
	unsigned int from_c = index_array[from_ed];
	unsigned int to_c = index_array[to_ed];
	parc = allocatePreArc ( to_c );
	parc->multiplicity = weight;
	parc->next = contig_array[from_c].arcs;
	contig_array[from_c].arcs = parc;
}

/*************************************************
 Function:
    loadArcs
 Description:
    Loads arc information of contigs and calculates the average weight of arcs.
 Input:
    1. graphfile:       prefix of graph file
 Output:
    None.
 Return:
    The average weight of arcs.
 *************************************************/
static unsigned int loadArcs ( char * graphfile )
{
	FILE * fp;
	char name[256], line[1024];
	unsigned int target, weight;
	unsigned int from_ed;
	char * seg;
	unsigned int avg_weight = 0, weight_sum = 0, arc_num = 0;
	sprintf ( name, "%s.Arc", graphfile );
	fp = ckopen ( name, "r" );
	createPreArcMemManager ();
	arcCounter = 0;

	while ( fgets ( line, sizeof ( line ), fp ) != NULL )
	{
		seg = strtok ( line, " " );
		from_ed = atoi ( seg );

		//printf("%d\n",from_ed);
		while ( ( seg = strtok ( NULL, " " ) ) != NULL )
		{
			target = atoi ( seg );
			seg = strtok ( NULL, " " );
			weight = atoi ( seg );
			add1Arc ( from_ed, target, weight );

			if ( !contig_array[index_array[from_ed]].mask && !contig_array[index_array[target]].mask )
			{
				weight_sum += weight;
				++arc_num;
			}
		}
	}

	if ( arc_num )
	{
		avg_weight = weight_sum / arc_num;
	}

	fprintf ( stderr, "%lld arc(s) loaded, average weight is %u.\n", arcCounter, avg_weight );
	fclose ( fp );
	return avg_weight;
}

/*************************************************
 Function:
    loadContig
 Description:
    Loads contigs sequence.
 Input:
    1. graphfile:       prefix of graph file
 Output:
    None.
 Return:
    None.
 *************************************************/
void loadContig ( char * graphfile )
{
	char c, name[256], line[1024], *tightSeq = NULL;
	FILE * fp;
	int n = 0, length, index = -1, edgeno;
	unsigned int i;
	unsigned int newIndex;
	sprintf ( name, "%s.contig", graphfile );
	fp = ckopen ( name, "r" );

	while ( fgets ( line, sizeof ( line ), fp ) != NULL )
	{
		if ( line[0] == '>' )
		{
			if ( index >= 0 )
			{
				newIndex = index_array[edgeno];
				contig_array[newIndex].seq = tightSeq;
			}

			n = 0;
			index++;
			sscanf ( line + 1, "%d %s %d", &edgeno, name, &length );
			//printf("contig %d, length %d\n",edgeno,length);
			tightSeq = ( char * ) ckalloc ( ( length / 4 + 1 ) * sizeof ( char ) );
		}
		else
		{
			for ( i = 0; i < strlen ( line ); i++ )
			{
				if ( line[i] >= 'a' && line[i] <= 'z' )
				{
					c = base2int ( line[i] - 'a' + 'A' );
					writeChar2tightString ( c, tightSeq, n++ );
				}
				else if ( line[i] >= 'A' && line[i] <= 'Z' )
				{
					c = base2int ( line[i] );
					writeChar2tightString ( c, tightSeq, n++ );
				}
			}
		}
	}

	if ( index >= 0 )
	{
		newIndex = index_array[edgeno];
		contig_array[newIndex].seq = tightSeq;
	}

	fprintf ( stderr, "%d contig(s) loaded.\n", index + 1 );
	fclose ( fp );
	//printf("the %dth contig with index 107\n",index);
}

void freeContig_array ()
{
	if ( !contig_array )
	{
		return;
	}

	unsigned int i;

	for ( i = 1; i <= num_ctg; i++ )
	{
		if ( contig_array[i].seq )
		{
			free ( ( void * ) contig_array[i].seq );
		}

		if ( contig_array[i].closeReads )
		{
			freeStack ( contig_array[i].closeReads );
		}
	}

	free ( ( void * ) contig_array );
	contig_array = NULL;
}

/*
void loadCvg(char *graphfile)
{
    char name[256],line[1024];
    FILE *fp;
    int cvg;
    unsigned int newIndex,edgeno,bal_ctg;

    sprintf(name,"%s.contigCVG",graphfile);
    fp = fopen(name,"r");
    if(!fp){
        printf("contig coverage file %s is not found!\n",name);
        return;
    }

    while(fgets(line,sizeof(line),fp)!=NULL){
        if(line[0]=='>'){
            sscanf(line+1,"%d %d",&edgeno,&cvg);
            newIndex = index_array[edgeno];
            cvg = cvg <= 255 ? cvg:255;
            contig_array[newIndex].multi = cvg;
            bal_ctg = getTwinCtg(newIndex);
            contig_array[bal_ctg].multi= cvg;
        }
    }
    fclose(fp);
}
*/

/*
 * loadPreGraph.c
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
#include "zlib.h"

static void loadPreArcs ( char * graphfile );

int cmp_vertex ( const void * a, const void * b )
{
	VERTEX * A, *B;
	A = ( VERTEX * ) a;
	B = ( VERTEX * ) b;

	if ( KmerLarger ( A->kmer, B->kmer ) )
	{
		return 1;
	}
	else if ( KmerEqual ( A->kmer, B->kmer ) )
	{
		return 0;
	}
	else
	{
		return -1;
	}
}

/*************************************************
Function:
    loadVertex
Description:
    1. Loads first kmer and the last kmer of the edge.
    2. Sorts kmers.
Input:
    1. graphfile:       the input prefix
Output:
    None.
Return:
    None.
*************************************************/
void loadVertex ( char * graphfile )
{
	char name[256], line[256];
	FILE * fp;
	Kmer word, bal_word, temp;
	int num_kmer, i;
	char ch;
	sprintf ( name, "%s.preGraphBasic", graphfile );
	fp = ckopen ( name, "r" );

	while ( fgets ( line, sizeof ( line ), fp ) != NULL )
	{
		if ( line[0] == 'V' )
		{
			sscanf ( line + 6, "%d %c %d", &num_kmer, &ch, &overlaplen );
			fprintf ( stderr, "There are %d kmer(s) in vertex file.\n", num_kmer );
		}
		else if ( line[0] == 'E' )
		{
			sscanf ( line + 5, "%d", &num_ed );
			fprintf ( stderr, "There are %d edge(s) in edge file.\n", num_ed );
		}
		else if ( line[0] == 'M' )
		{
			sscanf ( line, "MaxReadLen %d MinReadLen %d MaxNameLen %d", &maxReadLen, &minReadLen, &maxNameLen );
		}
		else if ( line[0] == 'B' )
		{
			if ( line[7] == 'V' )
			{
				sscanf ( line, "Backup VERTEX %d %c %d", &num_kmer, &ch, &overlaplen );
				fprintf ( stderr, "Backup there are %d kmer(s) in vertex file.\n", num_kmer );
			}
			else if ( line[7] == 'E' )
			{
				sscanf ( line, "Backup EDGEs %d", &num_ed );
				fprintf ( stderr, "Backup there are %d edge(s) in edge file.\n", num_ed );
			}
			else if ( line[7] == 'M' )
			{
				sscanf ( line, "Backup MaxReadLen %d MinReadLen %d MaxNameLen %d", &maxReadLen, &minReadLen, &maxNameLen );
			}
		}
	}

	fclose ( fp );
	vt_array = ( VERTEX * ) ckalloc ( ( 4 * num_kmer ) * sizeof ( VERTEX ) );
	num_kmer_limit = 4 * num_kmer;
	sprintf ( name, "%s.vertex", graphfile );
	fp = ckopen ( name, "r" );

	for ( i = 0; i < num_kmer; i++ )
	{
#ifdef MER127
		fscanf ( fp, "%llx %llx %llx %llx", & ( word.high1 ), & ( word.low1 ), & ( word.high2 ), & ( word.low2 ) );
#else
		fscanf ( fp, "%llx %llx", & ( word.high ), & ( word.low ) );
#endif
		bal_word = reverseComplement ( word, overlaplen );

		if ( KmerSmaller ( word, bal_word ) )
		{
			vt_array[i].kmer = word;
		}
		else
		{
			vt_array[i].kmer = bal_word;
		}
	}

	temp = vt_array[num_kmer - 1].kmer;
	qsort ( &vt_array[0], num_kmer, sizeof ( vt_array[0] ), cmp_vertex );
	fprintf ( stderr, "Kmers sorted.\n" );
	fclose ( fp );

	for ( i = 0; i < num_kmer; i++ )
	{
		bal_word = reverseComplement ( vt_array[i].kmer, overlaplen );
		vt_array[i + num_kmer].kmer = bal_word;
	}

	num_vt = num_kmer;
}

/*************************************************
Function:
    bisearch
Description:
    Uses binary search to get the index of a kmer in array.
Input:
    1. vts:     the array of the kmer
    2. num:     the count of the total kmers
    3. target:      the kmer to search
Output:
    None.
Return:
    The kmer's index in array.
*************************************************/
int bisearch ( VERTEX * vts, int num, Kmer target )
{
	int mid, low, high;
	low = 0;
	high = num - 1;

	while ( low <= high )
	{
		mid = ( low + high ) / 2;

		if ( KmerEqual ( vts[mid].kmer, target ) )
		{
			break;
		}
		else if ( KmerLarger ( target, vts[mid].kmer ) )
		{
			low = mid + 1;
		}
		else
		{
			high = mid - 1;
		}
	}

	if ( low <= high )
	{
		return mid;
	}
	else
	{
		return -1;
	}
}

/*************************************************
Function:
    kmer2vt
Description:
    Searchs the index of a kmer in array.
Input:
    1. kmer:        the kmer
Output:
    None.
Return:
    The index of the kmer in the array.
*************************************************/
int kmer2vt ( Kmer kmer )
{
	Kmer bal_word;
	int vt_id;
	bal_word = reverseComplement ( kmer, overlaplen );

	if ( KmerSmaller ( kmer, bal_word ) )
	{
		vt_id = bisearch ( &vt_array[0], num_vt, kmer );

		if ( vt_id < 0 )
		{
			fprintf ( stderr, "There is no vertex for kmer " );
			PrintKmer ( stderr, kmer );
			fprintf ( stderr, " .\n" );
			/*
			#ifdef MER127
			fprintf (stderr,"There is not the vertex for kmer %llx %llx %llx %llx.\n", kmer.high1, kmer.low1, kmer.high2, kmer.low2);
			#else
			fprintf (stderr,"There is not the vertex for kmer %llx %llx.\n", kmer.high, kmer.low);
			#endif
			*/
		}

		return vt_id;
	}
	else
	{
		vt_id = bisearch ( &vt_array[0], num_vt, bal_word );

		if ( vt_id >= 0 )
		{
			vt_id += num_vt;
		}
		else
		{
			fprintf ( stderr, "There is no vertex for kmer " );
			PrintKmer ( stderr, kmer );
			fprintf ( stderr, " .\n" );
			/*
			#ifdef MER127
			fprintf (stderr,"There is not the vertex for kmer %llx %llx %llx %llx.\n", kmer.high1, kmer.low1, kmer.high2, kmer.low2);
			#else
			fprintf (stderr,"There is not the vertex for kmer %llx %llx.\n", kmer.high, kmer.low);
			#endif
			*/
		}

		return vt_id;
	}
}

/*************************************************
Function:
    buildReverseComplementEdge
Description:
    Creates  reverse complementary edge of a edge.
Input:
    1. edgeno:      the edge index
Output:
    None.
Return:
    None.
*************************************************/
#ifdef MER127
static void buildReverseComplementEdge ( unsigned int edgeno )
{
	int length = edge_array[edgeno].length;
	int i, index = 0;
	char * sequence, ch, *tightSeq;
	Kmer kmer = vt_array[edge_array[edgeno].from_vt].kmer;
	sequence = ( char * ) ckalloc ( ( overlaplen + length ) * sizeof ( char ) );
	int bit1, bit2, bit3, bit4;

	if ( overlaplen < 32 )
	{
		bit4 = overlaplen;
		bit3 = 0;
		bit2 = 0;
		bit1 = 0;
	}

	if ( overlaplen >= 32 && overlaplen < 64 )
	{
		bit4 = 32;
		bit3 = overlaplen - 32;
		bit2 = 0;
		bit1 = 0;
	}

	if ( overlaplen >= 64 && overlaplen < 96 )
	{
		bit4 = 32;
		bit3 = 32;
		bit2 = overlaplen - 64;
		bit1 = 0;
	}

	if ( overlaplen >= 96 && overlaplen < 128 )
	{
		bit4 = 32;
		bit3 = 32;
		bit2 = 32;
		bit1 = overlaplen - 96;
	}

	for ( i = bit1 - 1; i >= 0; i-- )
	{
		ch = kmer.high1 & 0x3;
		kmer.high1 >>= 2;
		sequence[i] = ch;
	}

	for ( i = bit2 - 1; i >= 0; i-- )
	{
		ch = kmer.low1 & 0x3;
		kmer.low1 >>= 2;
		sequence[i + bit1] = ch;
	}

	for ( i = bit3 - 1; i >= 0; i-- )
	{
		ch = kmer.high2 & 0x3;
		kmer.high2 >>= 2;
		sequence[i + bit1 + bit2] = ch;
	}

	for ( i = bit4 - 1; i >= 0; i-- )
	{
		ch = kmer.low2 & 0x3;
		kmer.low2 >>= 2;
		sequence[i + bit1 + bit2 + bit3] = ch;
	}

	for ( i = 0; i < length; i++ )
	{
		sequence[i + overlaplen] = getCharInTightString ( edge_array[edgeno].seq, i );
	}

	tightSeq = ( char * ) ckalloc ( ( length / 4 + 1 ) * sizeof ( char ) );

	for ( i = length - 1; i >= 0; i-- )
	{
		writeChar2tightString ( int_comp ( sequence[i] ), tightSeq, index++ );
	}

	edge_array[edgeno + 1].length = length;
	edge_array[edgeno + 1].cvg = edge_array[edgeno].cvg;
	kmer = vt_array[edge_array[edgeno].from_vt].kmer;
	edge_array[edgeno + 1].to_vt = kmer2vt ( reverseComplement ( kmer, overlaplen ) );
	kmer = vt_array[edge_array[edgeno].to_vt].kmer;
	edge_array[edgeno + 1].from_vt = kmer2vt ( reverseComplement ( kmer, overlaplen ) );
	edge_array[edgeno + 1].seq = tightSeq;
	edge_array[edgeno + 1].bal_edge = 0;
	edge_array[edgeno + 1].rv = NULL;
	edge_array[edgeno + 1].arcs = NULL;
	edge_array[edgeno + 1].flag = 0;
	edge_array[edgeno + 1].deleted = 0;
	free ( ( void * ) sequence );
}
#else

/*************************************************
Function:
    buildReverseComplementEdge
Description:
    Creates reverse complementary edge of a edge.
Input:
    1. edgeno:      the edge index
Output:
    None.
Return:
    None.
*************************************************/
static void buildReverseComplementEdge ( unsigned int edgeno )
{
	int length = edge_array[edgeno].length;
	int i, index = 0;
	char * sequence, ch, *tightSeq;
	Kmer kmer = vt_array[edge_array[edgeno].from_vt].kmer;
	sequence = ( char * ) ckalloc ( ( overlaplen + length ) * sizeof ( char ) );
	int bit2 = overlaplen > 32 ? 32 : overlaplen;
	int bit1 = overlaplen > 32 ? overlaplen - 32 : 0;

	for ( i = bit1 - 1; i >= 0; i-- )
	{
		ch = kmer.high & 0x3;
		kmer.high >>= 2;
		sequence[i] = ch;
	}

	for ( i = bit2 - 1; i >= 0; i-- )
	{
		ch = kmer.low & 0x3;
		kmer.low >>= 2;
		sequence[i + bit1] = ch;
	}

	for ( i = 0; i < length; i++ )
	{
		sequence[i + overlaplen] = getCharInTightString ( edge_array[edgeno].seq, i );
	}

	tightSeq = ( char * ) ckalloc ( ( length / 4 + 1 ) * sizeof ( char ) );

	for ( i = length - 1; i >= 0; i-- )
	{
		writeChar2tightString ( int_comp ( sequence[i] ), tightSeq, index++ );
	}

	edge_array[edgeno + 1].length = length;
	edge_array[edgeno + 1].cvg = edge_array[edgeno].cvg;
	kmer = vt_array[edge_array[edgeno].from_vt].kmer;
	edge_array[edgeno + 1].to_vt = kmer2vt ( reverseComplement ( kmer, overlaplen ) );
	kmer = vt_array[edge_array[edgeno].to_vt].kmer;
	edge_array[edgeno + 1].from_vt = kmer2vt ( reverseComplement ( kmer, overlaplen ) );
	edge_array[edgeno + 1].seq = tightSeq;
	edge_array[edgeno + 1].bal_edge = 0;
	edge_array[edgeno + 1].rv = NULL;
	edge_array[edgeno + 1].arcs = NULL;
	edge_array[edgeno + 1].flag = 0;
	edge_array[edgeno + 1].deleted = 0;
	free ( ( void * ) sequence );
}
#endif

/*************************************************
Function:
    loadEdge
Description:
    1. Loads the info of the edges.
    2. Loads the links between the edges.
Input:
    1. graphfile:       the input prefix
Output:
    None.
Return:
    None.
*************************************************/
void loadEdge ( char * graphfile )
{
	char c, name[256], line[1024], str[32];
	char * tightSeq = NULL;
	gzFile * fp;
	Kmer from_kmer, to_kmer;
	int n = 0, i, length, cvg, index = -1, bal_ed, edgeno;
	int linelen;
	unsigned int j;
	sprintf ( name, "%s.edge.gz", graphfile );
	fp = gzopen ( name, "r" );
	num_ed_limit = 1.2 * num_ed;
	edge_array = ( EDGE * ) ckalloc ( ( num_ed_limit + 1 ) * sizeof ( EDGE ) );

	for ( j = num_ed + 1; j <= num_ed_limit; j++ )
	{
		edge_array[j].seq = NULL;
	}

	while ( gzgets ( fp, line, sizeof ( line ) ) != NULL )
	{
		if ( line[0] == '>' )
		{
			if ( index >= 0 )
			{
				edgeno = index + 1;
				edge_array[edgeno].length = length;
				edge_array[edgeno].cvg = cvg;
				edge_array[edgeno].from_vt = kmer2vt ( from_kmer );
				edge_array[edgeno].to_vt = kmer2vt ( to_kmer );
				edge_array[edgeno].seq = tightSeq;
				edge_array[edgeno].bal_edge = bal_ed + 1;
				edge_array[edgeno].rv = NULL;
				edge_array[edgeno].arcs = NULL;
				edge_array[edgeno].flag = 0;
				edge_array[edgeno].deleted = 0;

				if ( bal_ed )
				{
					buildReverseComplementEdge ( edgeno );
					index++;
				}
			}

			n = 0;
			index++;
#ifdef MER127
			sscanf ( line + 7, "%d,%llx %llx %llx %llx,%llx %llx %llx %llx,%s %d,%d",
			         &length, & ( from_kmer.high1 ), & ( from_kmer.low1 ), & ( from_kmer.high2 ), & ( from_kmer.low2 ), & ( to_kmer.high1 ), & ( to_kmer.low1 ),
			         & ( to_kmer.high2 ), & ( to_kmer.low2 ), str, &cvg, &bal_ed );
#else
			sscanf ( line + 7, "%d,%llx %llx,%llx %llx,%s %d,%d", &length, & ( from_kmer.high ), & ( from_kmer.low ), & ( to_kmer.high ), & ( to_kmer.low ), str, &cvg, &bal_ed );
#endif
			tightSeq = ( char * ) ckalloc ( ( length / 4 + 1 ) * sizeof ( char ) );
		}
		else
		{
			linelen = strlen ( line );

			for ( i = 0; i < linelen; i++ )
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
		edgeno = index + 1;
		edge_array[edgeno].length = length;
		edge_array[edgeno].cvg = cvg;
		edge_array[edgeno].from_vt = kmer2vt ( from_kmer );
		edge_array[edgeno].to_vt = kmer2vt ( to_kmer );
		edge_array[edgeno].seq = tightSeq;
		edge_array[edgeno].bal_edge = bal_ed + 1;

		if ( bal_ed )
		{
			buildReverseComplementEdge ( edgeno );
			index++;
		}
	}

	fprintf ( stderr, "%d edge(s) input.\n", index + 1 );
	gzclose ( fp );
	createArcMemo ();
	loadPreArcs ( graphfile );
}

unsigned int getTwinEdge ( unsigned int edgeno )
{
	return edgeno + edge_array[edgeno].bal_edge - 1;
}

boolean EdSmallerThanTwin ( unsigned int edgeno )
{
	return edge_array[edgeno].bal_edge > 1;
}

boolean EdLargerThanTwin ( unsigned int edgeno )
{
	return edge_array[edgeno].bal_edge < 1;
}

boolean EdSameAsTwin ( unsigned int edgeno )
{
	return edge_array[edgeno].bal_edge == 1;
}

/*************************************************
Function:
    add1Arc
Description:
    Updates the arc information between two edges.
Input:
    1. from_ed:     the first edge
    2. to_ed:           the second edge
    3. weight:      the weight of the arc to add
Output:
    None.
Return:
    None.
*************************************************/
static void add1Arc ( unsigned int from_ed, unsigned int to_ed, unsigned int weight )
{
	if ( edge_array[from_ed].to_vt != edge_array[to_ed].from_vt )
	{
		//fprintf(stderr,"add1Arc: inconsistant joins\n");
		return;
	}

	unsigned int bal_fe = getTwinEdge ( from_ed );
	unsigned int bal_te = getTwinEdge ( to_ed );

	if ( from_ed > num_ed || to_ed > num_ed || bal_fe > num_ed || bal_te > num_ed )
	{
		return;
	}

	ARC * parc, *bal_parc;
	//both arcs already exist
	parc = getArcBetween ( from_ed, to_ed );

	if ( parc )
	{
		bal_parc = parc->bal_arc;
		parc->multiplicity += weight;
		bal_parc->multiplicity += weight;
		return;
	}

	//create new arcs
	parc = allocateArc ( to_ed );
	parc->multiplicity = weight;
	parc->prev = NULL;

	if ( edge_array[from_ed].arcs )
	{
		edge_array[from_ed].arcs->prev = parc;
	}

	parc->next = edge_array[from_ed].arcs;
	edge_array[from_ed].arcs = parc;

	// A->A'
	if ( bal_te == from_ed )
	{
		//printf("preArc from A to A'\n");
		parc->bal_arc = parc;
		parc->multiplicity += weight;
		return;
	}

	bal_parc = allocateArc ( bal_fe );
	bal_parc->multiplicity = weight;
	bal_parc->prev = NULL;

	if ( edge_array[bal_te].arcs )
	{
		edge_array[bal_te].arcs->prev = bal_parc;
	}

	bal_parc->next = edge_array[bal_te].arcs;
	edge_array[bal_te].arcs = bal_parc;
	//link them to each other
	parc->bal_arc = bal_parc;
	bal_parc->bal_arc = parc;
}

/*************************************************
Function:
    loadPreArcs
Description:
    Loads the info of pre-arcs.
Input:
    1. graphfile:       the input prefix
Output:
    None.
Return:
    None.
*************************************************/
void loadPreArcs ( char * graphfile )
{
	FILE * fp;
	char name[256], line[1024];
	unsigned int target, weight;
	unsigned int from_ed;
	char * seg;
	sprintf ( name, "%s.preArc", graphfile );
	fp = ckopen ( name, "r" );
	arcCounter = 0;

	while ( fgets ( line, sizeof ( line ), fp ) != NULL )
	{
		seg = strtok ( line, " " );
		from_ed = atoi ( seg );

		while ( ( seg = strtok ( NULL, " " ) ) != NULL )
		{
			target = atoi ( seg );
			seg = strtok ( NULL, " " );
			weight = atoi ( seg );
			add1Arc ( from_ed, target, weight );
		}
	}

	fprintf ( stderr, "%lli pre-arcs loaded.\n", arcCounter );
	fclose ( fp );
}

void free_edge_array ( EDGE * ed_array, int ed_num )
{
	int i;

	for ( i = 1; i <= ed_num; i++ )
		if ( ed_array[i].seq )
		{
			free ( ( void * ) ed_array[i].seq );
		}

	free ( ( void * ) ed_array );
}

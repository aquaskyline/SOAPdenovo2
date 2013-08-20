/*
 * pregraph.c
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


static char shortrdsfile[256];   //the reads config file name ,see -s option
static char graphfile[256];      //the output prefix name ,see -o option
static int cutTips = 1;          //whether remove single tips or not. single tips , the tips starting from a kmer which coverage = 1

static void initenv ( int argc, char ** argv );
static void display_pregraph_usage ();


/*************************************************
Function:
    call_pregraph
Description:
    The main function for pregraph step . its processes are as below:
    1. Builds the kmer hash sets and remove the low coverage kmers.
    2. Removes the tips which length are no greater than 2*K.
    3. Builds edges by combining  linear kmers.
    4. Maps the reads back to edges and build preArcs (the connection between edges).
Input:
    @see display_pregraph_usage ()
Output:
    Below files:
    *.kmerFreq
    *.edge.gz
    *.vertex
    *.preArc
    *.preGraphBasic
    *.markOnEdge    (optional)
    *.path          (optional)
Return:
    Zero always
*************************************************/

int call_pregraph ( int argc, char ** argv )
{
	time_t start_t, stop_t, time_bef, time_aft;
	time ( &start_t );
	fprintf ( stderr, "\n********************\n" );
	fprintf ( stderr, "Pregraph\n" );
	fprintf ( stderr, "********************\n\n" );
	initenv ( argc, argv );

	if ( overlaplen % 2 == 0 )
	{
		overlaplen++;
		fprintf ( stderr, "K should be an odd number.\n" );
	}

	if ( overlaplen < 13 )
	{
		overlaplen = 13;
		fprintf ( stderr, "K should not be less than 13.\n" );
	}

#ifdef MER127
	else if ( overlaplen > 127 )
	{
		overlaplen = 127;
		fprintf ( stderr, "K should not be greater than 127.\n" );
	}

#else
	else if ( overlaplen > 63 )
	{
		overlaplen = 63;
		fprintf ( stderr, "K should not be greater than 63.\n" );
	}

#endif
	time ( &time_bef );
	prlRead2HashTable ( shortrdsfile, graphfile );
	time ( &time_aft );
	fprintf ( stderr, "Time spent on pre-graph construction: %ds.\n\n", ( int ) ( time_aft - time_bef ) );
	//  printf ("deLowKmer %d, deLowEdge %d\n", deLowKmer, deLowEdge);
	//  fprintf (stderr,"DeLowKmer %d\n", deLowKmer);

	//analyzeTips(hash_table, graphfile);
	if ( !deLowKmer && cutTips )
	{
		time ( &time_bef );
		removeSingleTips ();
		removeMinorTips ();
		time ( &time_aft );
		fprintf ( stderr, "Time spent on removing tips: %ds.\n\n", ( int ) ( time_aft - time_bef ) );
	}
	else
	{
		time ( &time_bef );
		removeMinorTips ();
		time ( &time_aft );
		fprintf ( stderr, "Time spent on removing tips: %ds.\n\n", ( int ) ( time_aft - time_bef ) );
	}

	initKmerSetSize = 0;
	//combine each linear part to an edge
	time ( &time_bef );
	kmer2edges ( graphfile );
	time ( &time_aft );
	fprintf ( stderr, "Time spent on constructing edges: %ds.\n\n", ( int ) ( time_aft - time_bef ) );
	//map read to edge one by one
	time ( &time_bef );
	prlRead2edge ( shortrdsfile, graphfile );
	time ( &time_aft );
	fprintf ( stderr, "Time spent on aligning reads: %ds.\n\n", ( int ) ( time_aft - time_bef ) );
	output_vertex ( graphfile );
	free_Sets ( KmerSets, thrd_num );
	free_Sets ( KmerSetsPatch, thrd_num );
	time ( &stop_t );
	fprintf ( stderr, "Overall time spent on constructing pre-graph: %dm.\n\n", ( int ) ( stop_t - start_t ) / 60 );
	return 0;
}


void initenv ( int argc, char ** argv )
{
	int copt;
	int inpseq, outseq;
	extern char * optarg;
	char temp[100];
	optind = 1;
	inpseq = outseq = 0;
	fprintf ( stderr, "Parameters: pregraph " );

	while ( ( copt = getopt ( argc, argv, "a:s:o:K:p:d:R" ) ) != EOF )
	{
		//printf("get option\n");
		switch ( copt )
		{
			case 's':
				fprintf ( stderr, "-s %s ", optarg );
				inpseq = 1;
				sscanf ( optarg, "%s", shortrdsfile );
				break;
			case 'o':
				fprintf ( stderr, "-o %s ", optarg );
				outseq = 1;
				sscanf ( optarg, "%s", graphfile );
				break;
			case 'K':
				fprintf ( stderr, "-K %s ", optarg );
				sscanf ( optarg, "%s", temp );
				overlaplen = atoi ( temp );
				break;
			case 'p':
				fprintf ( stderr, "-p %s ", optarg );
				sscanf ( optarg, "%s", temp );
				thrd_num = atoi ( temp );
				break;
			case 'R':
				repsTie = 1;
				fprintf ( stderr, "-R " );
				break;
			case 'd':
				fprintf ( stderr, "-d %s ", optarg );
				sscanf ( optarg, "%s", temp );
				deLowKmer = atoi ( temp ) >= 0 ? atoi ( temp ) : 0;
				break;
				/*
				            case 'D':
				                deLowEdge = 1;
				                break;
				*/
			case 'a':
				fprintf ( stderr, "-a %s ", optarg );
				initKmerSetSize = atoi ( optarg );
				break;
			default:

				if ( inpseq == 0 || outseq == 0 )
				{
					display_pregraph_usage ();
					exit ( -1 );
				}
		}
	}

	fprintf ( stderr, "\n\n" );

	if ( inpseq == 0 || outseq == 0 )
	{
		//printf("need more\n");
		display_pregraph_usage ();
		exit ( -1 );
	}
}

static void display_pregraph_usage ()
{
	fprintf ( stderr, "\npregraph -s configFile -o outputGraph [-R] [-K kmer -p n_cpu -a initMemoryAssumption -d KmerFreqCutoff]\n" );
	fprintf ( stderr, "  -s <string>      configFile: the config file of solexa reads\n" );
	fprintf ( stderr, "  -o <string>      outputGraph: prefix of output graph file name\n" );
#ifdef MER127
	fprintf ( stderr, "  -K <int>         kmer(min 13, max 127): kmer size, [23]\n" );
#else
	fprintf ( stderr, "  -K <int>         kmer(min 13, max 63): kmer size, [23]\n" );
#endif
	fprintf ( stderr, "  -p <int>         n_cpu: number of cpu for use, [8]\n" );
	fprintf ( stderr, "  -a <int>         initMemoryAssumption: memory assumption initialized to avoid further reallocation, unit GB, [0]\n" );
	fprintf ( stderr, "  -R (optional)    output extra information for resolving repeats in contig step, [NO]\n" );
	fprintf ( stderr, "  -d <int>         KmerFreqCutoff: kmers with frequency no larger than KmerFreqCutoff will be deleted, [0]\n" );
}

/*
 * contig.c
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
static void initenv ( int argc, char ** argv );
static void display_contig_usage ();
char shortrdsfile[256], graphfile[256];
static boolean repeatSolve;     //whether solve repeat or not
//static boolean keepReadFile = 0;  //whether keep tmp selected reads file or not
static boolean iter = 0;                //whether use multikmer or not
static boolean cleanBubble = 0;     //whether merge clean bubble before iterate
static int M = 1;           //merge bubble level
static int maxk = 0;        //max kmer of multikmer

/*************************************************
Function:
    call_heavygraph
Description:
    The main function for contig step . its processes are as below:
    1. Solve repeat
    2. Merge bubble(clean bubble is optional for multikmer)
    3. Remove weak edge and low coverage edge
    4. Cut tips
    5. Iterate multikmer(optional)
Input:
    @see display_contig_usage ()
Output:
    The files below:
    1. *.contig
    2. *.ContigIndex
    3. *.update.edge
    4. *.Arc
    5. *.read           [optional]
    6. *.preGraphBasic  [optional]
Return:
    None.
*************************************************/
int call_heavygraph ( int argc, char ** argv )
{
	time_t start_t, stop_t, time_bef, time_aft;
	time ( &start_t );
	boolean ret;
	fprintf ( stderr, "\n********************\n" );
	fprintf ( stderr, "Contig\n" );
	fprintf ( stderr, "********************\n\n" );
	initenv ( argc, argv );
	loadVertex ( graphfile );
	loadEdge ( graphfile );

	if ( repeatSolve )
	{
		ret = loadPathBin ( graphfile );
	}

	swapedge();
	sortedge();
	freshArc();

	if ( repeatSolve )
	{
		time ( &time_bef );

		//      ret = loadPathBin (graphfile);
		if ( ret )
		{
			solveReps ();
		}
		else
		{
			fprintf ( stderr, "Repeat solving can't be done...\n" );
		}

		time ( &time_aft );
		fprintf ( stderr, "Time spent on solving repeat: %ds.\n", ( int ) ( time_aft - time_bef ) );
	}

	//edgecvg_bar(edge_array,num_ed,graphfile,100);

	if ( !iter && M > 0 )
	{
		time ( &time_bef );
		bubblePinch ( 0.90, graphfile, M, 0, 1 );
		time ( &time_aft );
		fprintf ( stderr, "Time spent on pinching bubbles: %ds.\n", ( int ) ( time_aft - time_bef ) );
	}

	if ( iter && cleanBubble && M > 0 )
	{
		time ( &time_bef );
		clean = 1;
		long long oldpinCounter = 0;
		long long min = 10;
		int times = 0;

		while ( min >= 10 )
		{
			times++;

			if ( times >= 4 ) { break; }

			bubblePinch ( 0.90, graphfile, M, 1, 0 );
			min = pinCounter;
			fprintf ( stderr, "%lld clean bubbles merged.\n", pinCounter );
		}

		time ( &time_aft );
		fprintf ( stderr, "Time spent on pinching clean bubbles: %ds.\n", ( int ) ( time_aft - time_bef ) );
		clean = 0;
	}

	if ( deLowEdge )
	{
		removeWeakEdges ( 2 * overlaplen, 1 );
		removeLowCovEdges ( 2 * overlaplen, deLowEdge, !iter );
	}

	cutTipsInGraph ( 0, 0, !iter );

	if ( iter )
	{
		Iterate ( shortrdsfile, graphfile, maxk, M ); //keepReadFile,

		if ( M > 0 )
		{
			time ( &time_bef );
			bubblePinch ( 0.90, graphfile, M, 1, 0 );
			time ( &time_aft );
			fprintf ( stderr, "Time spent on pinching bubbles: %ds.\n", ( int ) ( time_aft - time_bef ) );
		}

		freshpreGraphBasic ( iter, maxk, graphfile );
	}

	//output_graph(graphfile);
	output_contig ( edge_array, num_ed, graphfile, overlaplen + 1 );
	output_updated_edges ( graphfile );
	output_heavyArcs ( graphfile );

	if ( vt_array )
	{
		free ( ( void * ) vt_array );
		vt_array = NULL;
	}

	if ( edge_array )
	{
		free_edge_array ( edge_array, num_ed_limit );
		edge_array = NULL;
	}

	destroyArcMem ();
	time ( &stop_t );
	fprintf ( stderr, "\nTime spent on constructing contig: %dm.\n\n", ( int ) ( stop_t - start_t ) / 60 );
	return 0;
}


/*****************************************************************************
 * Parse command line switches
 *****************************************************************************/
void initenv ( int argc, char ** argv )
{
	int copt;
	int inpseq, outseq;
	extern char * optarg;
	char temp[100];
	inpseq = outseq = repeatSolve = iter = cleanBubble = 0;//keepReadFile =
	optind = 1;
	fprintf ( stderr, "Parameters: contig " );

	while ( ( copt = getopt ( argc, argv, "g:M:D:Rs:m:p:e:E" ) ) != EOF ) // r
	{
		switch ( copt )
		{
			case 'M':
				fprintf ( stderr, "-M %s ", optarg );
				sscanf ( optarg, "%s", temp );
				M = atoi ( temp );
				break;
			case 'D':
				fprintf ( stderr, "-D %s ", optarg );
				sscanf ( optarg, "%s", temp );
				deLowEdge = atoi ( temp ) >= 0 ? atoi ( temp ) : 0;
				break;
			case 'g':
				fprintf ( stderr, "-g %s ", optarg );
				inGraph = 1;
				sscanf ( optarg, "%s", graphfile );
				break;
			case 'R':
				repeatSolve = 1;
				fprintf ( stderr, "-R " );
				break;
			case 's':
				fprintf ( stderr, "-s %s ", optarg );
				inpseq = 1;
				sscanf ( optarg, "%s", shortrdsfile );
				break;
			case 'm':
				fprintf ( stderr, "-m %s ", optarg );
				iter = 1;
				sscanf ( optarg, "%s", temp );
				maxk = atoi ( temp );
				break;
				/*
				case 'r':
				    keepReadFile = 1;
				    fprintf(stderr, "-r ");
				    break;
				    */
			case 'e':
				fprintf ( stderr, "-e %s ", optarg );
				sscanf ( optarg, "%s", temp );
				arcfilter = atoi ( temp );
				break;
			case 'p':
				fprintf ( stderr, "-p %s ", optarg );
				sscanf ( optarg, "%s", temp );
				thrd_num = atoi ( temp );
				break;
			case 'E':
				cleanBubble = 1;
				fprintf ( stderr, "-E " );
				break;
			default:

				if ( ( iter && inpseq == 0 ) || inGraph == 0 )
				{
					display_contig_usage ();
					exit ( -1 );
				}
		}
	}

	fprintf ( stderr, "\n\n" );

	if ( iter )
	{
		if ( maxk % 2 == 0 )
		{
			maxk++;
			fprintf ( stderr, "Max K should be an odd number, change to %d.\n", maxk );
		}

		if ( maxk < 13 )
		{
			maxk = 13;
			fprintf ( stderr, "Max K should not be less than 13, change to %d.\n", maxk );
		}

#ifdef MER127
		else if ( maxk > 127 )
		{
			maxk = 127;
			fprintf ( stderr, "Max K should not be greater than 127, change to %d.\n", maxk );
		}

#else
		else if ( maxk > 63 )
		{
			maxk = 63;
			fprintf ( stderr, "Max K should not be greater than 63, change to %d.\n", maxk );
		}

#endif

		if ( maxk <= overlaplen )
		{
			fprintf ( stderr, "Max K %d is not greater than overlaplen %d.\n", maxk, overlaplen );
			display_contig_usage ();
			exit ( -1 );
		}
	}

	if ( ( iter && inpseq == 0 ) || inGraph == 0 )
	{
		display_contig_usage ();
		exit ( -1 );
	}
}

static void display_contig_usage ()
{
	fprintf ( stderr, "\ncontig -g InputGraph [-R] [-M mergeLevel -D EdgeCovCutoff] [-s readsInfoFile -m maxkmer -p n_cpu -r]\n" );
	fprintf ( stderr, "  -g <string>      inputGraph: prefix of input graph file names\n" );
	fprintf ( stderr, "  -R (optional)    resolve repeats using information generated in pregraph step, works only if -R is set in pregraph step too, [NO]\n" );
	fprintf ( stderr, "  -M <int>         mergeLevel(min 0, max 3): the strength of merging similar sequences during contiging, [1]\n" );
	fprintf ( stderr, "  -D <int>         EdgeCovCutoff: edges shorter than (2*K+1) with coverage no larger than EdgeCovCutoff will be deleted, [1]\n" );
	fprintf ( stderr, "  -e <int>         arcWeight: two edges, between which the arc's weight is larger than arcWeight, will be linerized, [0]\n" );
	fprintf ( stderr, "  -m <int>         max k when using multi-kmer, and the parameters below are used along with multi-kmer, [NO]\n" );
	fprintf ( stderr, "      -s <string>      readsInfoFile:The file contains information of solexa reads(It's necessary when using multi-kmer)\n" );
	fprintf ( stderr, "      -p <int>         number of cpu, [8]\n" );
	fprintf ( stderr, "      -E (optional)    merge clean bubble before iterate, works only if -M is set when using multi-kmer, [NO]\n" );
	//  fprintf (stderr,"  -r (optional)    keep available read(*.read)\n");
}


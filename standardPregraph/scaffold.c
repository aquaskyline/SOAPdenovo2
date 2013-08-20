/*
 * scaffold.c
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
static void display_scaff_usage ();

static boolean LINK, SCAFF;
static char graphfile[256];

/*************************************************
Function:
    call_scaffold
Description:
    This is the main function of 'scaff' step. It includes the following steps:
    1. Reads required input information.
    2. Constructs scaffolds.
    3. Fills gaps if required, then output result.
Input:
    See @display_scaff_usage ().
Output:
    Scaffolds results and related information.
    1. *.scafSeq
    2. *.scafStatistics
    3. *.contigPosInscaff
    4. *.bubbleInScaff
    5. *.gapSeq
    6. *.links
    7. *.scaf
    8. *.scaf_gap
    9. *.newContigIndex
Return:
    0 if exit normally.
*************************************************/
int call_scaffold ( int argc, char ** argv )
{
	time_t start_t, stop_t, time_bef, time_aft;
	fprintf ( stderr, "\n********************\n" );
	fprintf ( stderr, "Scaff\n" );
	fprintf ( stderr, "********************\n\n" );
	ctg_short = 0;
	time ( &start_t );
	initenv ( argc, argv );
	checkFiles4Scaff ( graphfile );
	loadPEgrads ( graphfile );
	time ( &time_bef );
	loadUpdatedEdges ( graphfile );
	time ( &time_aft );
	fprintf ( stderr, "Time spent on loading updated edges: %ds.\n\n", ( int ) ( time_aft - time_bef ) );

	if ( !SCAFF )
	{
		time ( &time_bef );
		PE2Links ( graphfile );
		time ( &time_aft );
		fprintf ( stderr, "Time spent on loading paired-end reads information: %ds.\n", ( int ) ( time_aft - time_bef ) );
		time ( &time_bef );
		fprintf ( stderr, "\n*****************************************************\nStart to construct scaffolds.\n" );
		Links2Scaf ( graphfile );
		time ( &time_aft );
		fprintf ( stderr, "Time spent on constructing scaffolds: %ds.\n", ( int ) ( time_aft - time_bef ) );
		scaffolding ( 100, graphfile );
	}

	prlReadsCloseGap ( graphfile );
	ScafStat ( 100, graphfile );
	free_pe_mem ();

	if ( index_array )
	{
		free ( ( void * ) index_array );
	}

	freeContig_array ();
	destroyPreArcMem ();
	destroyConnectMem ();
	deleteCntLookupTable ();
	time ( &stop_t );
	fprintf ( stderr, "\nOverall time spent on constructing scaffolds: %dm.\n", ( int ) ( stop_t - start_t ) / 60 );
	return 0;
}

/*****************************************************************************
 * Parse command line switches
 *****************************************************************************/

void initenv ( int argc, char ** argv )
{
	int copt;
	int inpseq;
	extern char * optarg;
	char temp[256];
	inpseq = 0;
	LINK = 0;
	SCAFF = 0;
	optind = 1;
	fprintf ( stderr, "Parameters: scaff " );

	while ( ( copt = getopt ( argc, argv, "g:L:p:G:N:c:C:b:B:FzuSVw" ) ) != EOF )
	{
		switch ( copt )
		{
			case 'g':
				fprintf ( stderr, "-g %s ", optarg );
				inGraph = 1;
				sscanf ( optarg, "%s", graphfile );
				break;
			case 'G':
				fprintf ( stderr, "-G %s ", optarg );
				sscanf ( optarg, "%s", temp );
				GLDiff = atoi ( temp );
				break;
			case 'L':
				fprintf ( stderr, "-L %s ", optarg );
				sscanf ( optarg, "%s", temp );
				ctg_short = atoi ( temp );
				break;
			case 'N':
				fprintf ( stderr, "-N %s ", optarg );
				sscanf ( optarg, "%s", temp );
				known_genome_size = atoi ( temp );
				break;
			case 'F':
				fillGap = 1;
				fprintf ( stderr, "-F " );
				break;
			case 'u':
				maskRep = 0;
				fprintf ( stderr, "-u " );
				break;
			case 'S':
				SCAFF = 1;
				fprintf ( stderr, "-S " );
				break;
			case 'V':
				visual = 1;
				fprintf ( stderr, "-V " );
				break;
			case 'w':
				score_mask = 0;
				fprintf ( stderr, "-w " );
				break;
			case 'p':
				fprintf ( stderr, "-p %s ", optarg );
				sscanf ( optarg, "%s", temp );
				thrd_num = atoi ( temp );
				break;
			case 'c':
				fprintf ( stderr, "-c %s ", optarg );
				sscanf ( optarg, "%s", temp );
				cvg_low = atof ( temp ) > 0 ? atof ( temp ) : 0.0;
				break;
			case 'C':
				fprintf ( stderr, "-C %s ", optarg );
				sscanf ( optarg, "%s", temp );
				cvg_high = atof ( temp ) > 0 ? atof ( temp ) : 0.0;
				break;
			case 'b':
				fprintf ( stderr, "-b %s ", optarg );
				sscanf ( optarg, "%s", temp );
				ins_var_idx = atof ( temp ) > 1.0 ? atof ( temp ) : 0.0;
				break;
			case 'B':
				fprintf ( stderr, "-B %s ", optarg );
				sscanf ( optarg, "%s", temp );
				cvg4SNP = atof ( temp ) > 0 ? atof ( temp ) : 0.0;
				break;
			case 'z':
				COMPATIBLE_MODE = 1;
				fprintf ( stderr, "-m " );
				break;
			default:

				if ( inGraph == 0 )
				{
					display_scaff_usage ();
					exit ( -1 );
				}
		}
	}

	fprintf ( stderr, "\n\n" );

	if ( inGraph == 0 )
	{
		display_scaff_usage ();
		exit ( -1 );
	}
}

static void display_scaff_usage ()
{
	fprintf ( stderr, "\nscaff -g inputGraph [-F -z -u -S -w] [-G gapLenDiff -L minContigLen -c minContigCvg -C maxContigCvg -b insertSizeUpperBound -B bubbleCoverage -N genomeSize -p n_cpu]\n" );
	fprintf ( stderr, "  -g <string>        inputGraph: prefix of input graph file names\n" );
	fprintf ( stderr, "  -F (optional)      fill gaps in scaffold, [No]\n" );
	fprintf ( stderr, "  -z (optional)      use compatible mode to build scaffold with contig produced by Version 1.05, [No]\n" );
	fprintf ( stderr, "  -u (optional)      un-mask contigs with high/low coverage before scaffolding, [mask]\n" );
	fprintf ( stderr, "  -S (optional)      if scaffold structure exists, do gapfilling only(-F), [NO]\n" );
	fprintf ( stderr, "  -w (optional)      keep contigs weakly connected to other contigs in scaffold, [NO]\n" );
	fprintf ( stderr, "  -V (optional)      output information for Hawkeye to visualize the assembly, [NO]\n" );
	fprintf ( stderr, "  -G <int>           gapLenDiff: allowed length difference between estimated and filled gap, [50]\n" );
	fprintf ( stderr, "  -L <int>           minContigLen: shortest contig for scaffolding, [K+2]\n" );
	fprintf ( stderr, "  -c <float>         minContigCvg: minimum contig coverage (c*avgCvg), contigs shorter than 100bp with coverage smaller than c*avgCvg will be masked before scaffolding unless -u is set, [0.1]\n" );
	fprintf ( stderr, "  -C <float>         maxContigCvg: maximum contig coverage (C*avgCvg), contigs with coverage larger than C*avgCvg or contigs shorter than 100bp with coverage larger than 0.8*C*avgCvg will be masked before scaffolding unless -u is set, [2]\n" );
	fprintf ( stderr, "  -b <float>         insertSizeUpperBound: (b*avg_ins) will be used as upper bound of insert size for large insert size ( > 1000) when handling pair-end connections between contigs if b is set to larger than 1, [1.5]\n" );
	fprintf ( stderr, "  -B <float>         bubbleCoverage: remove contig with lower cvoerage in bubble structure if both contigs' coverage are smaller than bubbleCoverage*avgCvg, [0.6]\n" );
	fprintf ( stderr, "  -N <int>           genomeSize: genome size for statistics, [0]\n" );
	fprintf ( stderr, "  -p <int>           n_cpu: number of cpu for use, [8]\n" );
}

/*
 * map.c
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
static char shortrdsfile[256];
static char graphfile[256];

static void display_map_usage ();

/*************************************************
Function:
    getMinOverlap
Description:
    Parses file *.preGraphBasic, and get vertex number,
    K size, MaxReadLen, MinReadLen, MaxNameLen.
Input:
    1. gfile        the graph file prefix
Output:
    None.
Return:
    The kmer size.
*************************************************/
static int getMinOverlap ( char * gfile )
{
	char name[256], ch;
	FILE * fp;
	int num_kmer, overlaplen = 23;
	char line[1024];
	sprintf ( name, "%s.preGraphBasic", gfile );
	fp = fopen ( name, "r" );

	if ( !fp )
	{
		return overlaplen;
	}

	while ( fgets ( line, sizeof ( line ), fp ) != NULL )
	{
		if ( line[0] == 'V' )
		{
			sscanf ( line + 6, "%d %c %d", &num_kmer, &ch, &overlaplen );
		}
		else if ( line[0] == 'M' )
		{
			sscanf ( line, "MaxReadLen %d MinReadLen %d MaxNameLen %d", &maxReadLen, &minReadLen, &maxNameLen );
		}
	}

	fclose ( fp );
	return overlaplen;
}

/*************************************************
Function:
    call_align
Description:
    This is the main function for map step. It includes the following steps:
    1.Chops contig (>=k+2) into kmer sets, marks the kmer occured twice or more as deleted,
       records the unique kmers related contig's id and the position of the kmer on the contig.
    2.Maps reads to edges.
Input:
    @see display_map_usage()
Output:
    Below files:
    1. *.readOnCtg.gz
    2. *.readInGap.gz
Return:
    0 if exits normally.
*************************************************/

int call_align ( int argc, char ** argv )
{
	time_t start_t, stop_t, time_bef, time_aft;
	time ( &start_t );
	fprintf ( stderr, "\n********************\n" );
	fprintf ( stderr, "Map\n" );
	fprintf ( stderr, "********************\n\n" );
	initenv ( argc, argv );
	overlaplen = getMinOverlap ( graphfile );
#ifdef MER127

	if ( smallKmer > 12 && smallKmer < 128 && smallKmer % 2 == 1 )
	{
		deltaKmer = overlaplen - smallKmer;
		overlaplen = smallKmer;
	}

#else

	if ( smallKmer > 12 && smallKmer < 64 && smallKmer % 2 == 1 )
	{
		deltaKmer = overlaplen - smallKmer;
		overlaplen = smallKmer;
	}

#endif
	fprintf ( stderr, "Kmer size: %d.\n", overlaplen );
	time ( &time_bef );
	ctg_short = overlaplen + 2;
	fprintf ( stderr, "Contig length cutoff: %d.\n", ctg_short );
	prlContig2nodes ( graphfile, ctg_short );
	time ( &time_aft );
	fprintf ( stderr, "Time spent on graph construction: %ds.\n\n", ( int ) ( time_aft - time_bef ) );
	//map long read (asm_flags=4) to edge one by one
	time ( &time_bef );
	prlLongRead2Ctg ( shortrdsfile, graphfile );
	time ( &time_aft );
	fprintf ( stderr, "Time spent on aligning long reads: %ds.\n\n", ( int ) ( time_aft - time_bef ) );
	//map read to edge one by one
	time ( &time_bef );
	prlRead2Ctg ( shortrdsfile, graphfile );
	time ( &time_aft );
	fprintf ( stderr, "Time spent on aligning reads: %ds.\n\n", ( int ) ( time_aft - time_bef ) );
	free_Sets ( KmerSets, thrd_num );
	time ( &stop_t );
	fprintf ( stderr, "Overall time spent on alignment: %dm.\n\n", ( int ) ( stop_t - start_t ) / 60 );
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
	optind = 1;
	inpseq = outseq = 0;
	fprintf ( stderr, "Parameters: map " );

	while ( ( copt = getopt ( argc, argv, "s:g:K:p:k:f" ) ) != EOF )
	{
		//printf("get option\n");
		switch ( copt )
		{
			case 's':
				fprintf ( stderr, "-s %s ", optarg );
				inpseq = 1;
				sscanf ( optarg, "%s", shortrdsfile );
				break;
			case 'g':
				fprintf ( stderr, "-g %s ", optarg );
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
			case 'k':
				fprintf ( stderr, "-k %s ", optarg );
				sscanf ( optarg, "%s", temp );
				smallKmer = atoi ( temp );
				break;
			case 'f':
				fill = 1;
				fprintf ( stderr, "-f " );
				break;
			default:

				if ( inpseq == 0 || outseq == 0 )
				{
					display_map_usage ();
					exit ( 1 );
				}
		}
	}

	fprintf ( stderr, "\n\n" );

	if ( inpseq == 0 || outseq == 0 )
	{
		display_map_usage ();
		exit ( 1 );
	}
}

static void display_map_usage ()
{
	fprintf ( stderr, "\nmap -s configFile -g inputGraph [-f] [-p n_cpu -k kmer_R2C]\n" );
	fprintf ( stderr, "  -s <string>        configFile: the config file of solexa reads\n" );
	fprintf ( stderr, "  -g <string>        inputGraph: prefix of input graph file names\n" );
	fprintf ( stderr, "  -f (optional)      output gap related reads in map step for using SRkgf to fill gap, [NO]\n" );
	fprintf ( stderr, "  -p <int>           n_cpu: number of cpu for use, [8]\n" );
#ifdef MER127
	fprintf ( stderr, "  -k <int>           kmer_R2C(min 13, max 127): kmer size used for mapping read to contig, [K]\n" );
#else
	fprintf ( stderr, "  -k <int>           kmer_R2C(min 13, max 63): kmer size used for mapping read to contig, [K]\n" );
#endif
}

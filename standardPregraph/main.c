/*
 * main.c
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
#include "global.h"

extern int call_pregraph ( int arc, char ** argv );
extern int call_pregraph_sparse(int arc, char ** argv);
extern int call_heavygraph ( int arc, char ** argv );
extern int call_map2contig ( int arc, char ** argv );
extern int call_scaffold ( int arc, char ** argv );
extern int call_align ( int arc, char ** argv );

static void display_usage ();
static void display_all_usage ();
static void pipeline ( int argc, char ** argv );

/*************************************************
Function:
    main
Description:
    The main function. It includes four modules:
    1.pregraph
    2.contig
    3.map
    4.scaffold
Input:
    @see display_all_usage ()
Output:
    Below files:
    1. *.contig
    2. *.scafSeq
    etc.
Return:
    None.
*************************************************/
int main ( int argc, char ** argv )
{
	crc32c_Init();
	fprintf ( stderr, "\nVersion 2.04: released on July 13th, 2012\nCompile %s\t%s\n", __DATE__, __TIME__ );
	argc--;
	argv++;

	if ( argc == 0 )
	{
		display_usage ();
		return 0;
	}

	if ( strcmp ( "pregraph", argv[0] ) == 0 )
	{
		call_pregraph ( argc, argv );
	}
	else if(strcmp ( "sparse_pregraph", argv[0] ) == 0 ){
		call_pregraph_sparse ( argc, argv );

	}
	else if ( strcmp ( "contig", argv[0] ) == 0 )
	{
		call_heavygraph ( argc, argv );
	}
	else if ( strcmp ( "map", argv[0] ) == 0 )
	{
		call_align ( argc, argv );
	}
	//call_map2contig(argc,argv);
	else if ( strcmp ( "scaff", argv[0] ) == 0 )
	{
		call_scaffold ( argc, argv );
	}
	else if ( strcmp ( "all", argv[0] ) == 0 )
	{
		pipeline ( argc, argv );
	}
	else
	{
		display_usage ();
	}

	return 0;
}

static void display_usage ()
{
	fprintf ( stderr, "\nUsage: SOAPdenovo <command> [option]\n" );
	fprintf ( stderr, "    pregraph        construct kmer-graph\n" );
	fprintf ( stderr, "    sparse_pregraph construct sparse kmer-graph\n");
	fprintf ( stderr, "    contig          eliminate errors and output contigs\n" );
	fprintf ( stderr, "    map             map reads to contigs\n" );
	fprintf ( stderr, "    scaff           construct scaffolds\n" );
	fprintf ( stderr, "    all             do pregraph-contig-map-scaff in turn\n" );
}

static void pipeline ( int argc, char ** argv )
{
	char * options[32];
	unsigned char getK, getRfile, getOfile, getD, getDD, getL, getR, getP, getF, getf, getk, getu, getG, getc, getC, getb, getB, getN, getw, getV;
	unsigned char getm, getE; //getr,
	char readfile[256], outfile[256];
	char temp[128];
	char * name;
	int kmer = 0, cutoff_len = 0, ncpu = 0, lowK = 0, lowC = 0, kmer_small = 0, gap_diff = 0, genome_size = 0;
	float min_cvg = 0.0, max_cvg = 0.0, insert_size_bound = 0.0, bubble_coverage = 0.0;
	char kmer_s[16], len_s[16], ncpu_s[16], M_s[16], lowK_s[16], lowC_s[16], kmer_small_s[16], gap_diff_s[16], min_cvg_s[16], max_cvg_s[16], insert_size_bound_s[16], bubble_coverage_s[16], genome_size_s[16];
	int i, copt, index, M = 1;
	int maxk;
	char maxk_s[16];
	char arcfilter_s[16];
	extern char * optarg;
	time_t start_t, stop_t;
	time ( &start_t );
	getK = getRfile = getOfile = getD = getDD = getL = getR = getP = getF = getf = getk = getu = getG =  getc = getC = getb = getB = getN = getw = getm = getE = getV = 0;

	while ( ( copt = getopt ( argc, argv, "a:s:o:K:M:L:p:G:d:D:RuFk:fc:C:b:B:N:wm:e:EV" ) ) != EOF ) //r
	{
		switch ( copt )
		{
			case 's':
				getRfile = 1;
				sscanf ( optarg, "%s", readfile );
				break;
			case 'o':
				getOfile = 1;
				sscanf ( optarg, "%s", outfile );
				break;
			case 'K':
				getK = 1;
				sscanf ( optarg, "%s", temp );
				kmer = atoi ( temp );
				break;
			case 'G':
				getG = 1;
				sscanf ( optarg, "%s", temp );
				gap_diff = atoi ( temp );
				break;
			case 'M':
				sscanf ( optarg, "%s", temp );
				M = atoi ( temp );
				break;
			case 'p':
				getP = 1;
				sscanf ( optarg, "%s", temp );
				ncpu = atoi ( temp );
				break;
			case 'L':
				getL = 1;
				sscanf ( optarg, "%s", temp );
				cutoff_len = atoi ( temp );
				break;
			case 'R':
				getR = 1;
				break;
			case 'u':
				getu = 1;
				maskRep = 0;
				break;
			case 'd':
				getD = 1;
				sscanf ( optarg, "%s", temp );
				lowK = atoi ( temp );
				break;
			case 'D':
				getDD = 1;
				sscanf ( optarg, "%s", temp );
				lowC = atoi ( temp );
				break;
			case 'a':
				initKmerSetSize = atoi ( optarg );
				break;
			case 'F':
				getF = 1;
				break;
			case 'k':
				getk = 1;
				sscanf ( optarg, "%s", temp );
				kmer_small = atoi ( temp );
				break;
			case 'f':
				getf = 1;
				break;
			case 'c':
				getc = 1;
				sscanf ( optarg, "%s", temp );
				min_cvg = atof ( temp );
				break;
			case 'C':
				getC = 1;
				sscanf ( optarg, "%s", temp );
				max_cvg = atof ( temp );
				break;
			case 'b':
				getb = 1;
				sscanf ( optarg, "%s", temp );
				insert_size_bound = atof ( temp );
				break;
			case 'B':
				getB = 1;
				sscanf ( optarg, "%s", temp );
				bubble_coverage = atof ( temp );
				break;
			case 'N':
				getN = 1;
				sscanf ( optarg, "%s", temp );
				genome_size = atoi ( temp );
				break;
			case 'w':
				getw = 1;
				break;
			case 'm':
				getm = 1;
				sscanf ( optarg, "%s", temp );
				maxk = atoi ( temp );
				break;
				/*
				case 'r':
				getr = 1;
				break;
				*/
			case 'e':
				sscanf ( optarg, "%s", temp );
				arcfilter = atoi ( temp );
				break;
			case 'E':
				getE = 1;
				break;
			case 'V':
				getV = 1;
				break;
			default:

				if ( getRfile == 0 || getOfile == 0 )
				{
					display_all_usage ();
					exit ( -1 );
				}
		}
	}

	if ( getRfile == 0 || getOfile == 0 )
	{
		display_all_usage ();
		exit ( -1 );
	}

	if ( thrd_num < 1 )
	{
		thrd_num = 1;
	}

	// getK = getRfile = getOfile = getD = getL = getR = 0;
	name = "pregraph";
	index = 0;
	options[index++] = name;
	options[index++] = "-s";
	options[index++] = readfile;

	if ( getK )
	{
		options[index++] = "-K";
		sprintf ( kmer_s, "%d", kmer );
		options[index++] = kmer_s;
	}

	if ( getP )
	{
		options[index++] = "-p";
		sprintf ( ncpu_s, "%d", ncpu );
		options[index++] = ncpu_s;
	}

	if ( getD )
	{
		options[index++] = "-d";
		sprintf ( lowK_s, "%d", lowK );
		options[index++] = lowK_s;
	}

	if ( getR )
	{
		options[index++] = "-R";
	}

	options[index++] = "-o";
	options[index++] = outfile;
	/*
	for (i = 0; i < index; i++)
	{
	    fprintf (stderr,"%s ", options[i]);
	}

	fprintf (stderr,"\n");
	*/
	call_pregraph ( index, options );
	name = "contig";
	index = 0;
	options[index++] = name;
	options[index++] = "-g";
	options[index++] = outfile;
	options[index++] = "-M";
	sprintf ( M_s, "%d", M );
	options[index++] = M_s;

	if ( getR )
	{
		options[index++] = "-R";
	}

	if ( getDD )
	{
		options[index++] = "-D";
		sprintf ( lowC_s, "%d", lowC );
		options[index++] = lowC_s;
	}

	if ( getRfile )
	{
		options[index++] = "-s";
		options[index++] = readfile;
	}

	if ( getP )
	{
		options[index++] = "-p";
		sprintf ( ncpu_s, "%d", ncpu );
		options[index++] = ncpu_s;
	}

	if ( getm )
	{
		options[index++] = "-m";
		sprintf ( maxk_s, "%d", maxk );
		options[index++] = maxk_s;
	}

	/*
	if(getr){
	    options[index++] = "-r";
	}
	*/
	if ( getE )
	{
		options[index++] = "-E";
	}

	if ( arcfilter )
	{
		options[index++] = "-e";
		sprintf ( arcfilter_s, "%d", arcfilter );
		options[index++] = arcfilter_s;
	}

	/*
	for (i = 0; i < index; i++)
	{
	    fprintf (stderr,"%s ", options[i]);
	}

	fprintf (stderr,"\n");
	*/
	call_heavygraph ( index, options );
	name = "map";
	index = 0;
	options[index++] = name;
	options[index++] = "-s";
	options[index++] = readfile;
	options[index++] = "-g";
	options[index++] = outfile;

	if ( getP )
	{
		options[index++] = "-p";
		sprintf ( ncpu_s, "%d", ncpu );
		options[index++] = ncpu_s;
	}

	if ( getK )
	{
		options[index++] = "-K";
		sprintf ( kmer_s, "%d", kmer );
		options[index++] = kmer_s;
	}

	if ( getk )
	{
		options[index++] = "-k";
		sprintf ( kmer_small_s, "%d", kmer_small );
		options[index++] = kmer_small_s;
	}

	if ( getf )
	{
		options[index++] = "-f";
	}

	/*
	for (i = 0; i < index; i++)
	{
	    fprintf (stderr,"%s ", options[i]);
	}

	fprintf (stderr,"\n");
	*/
	call_align ( index, options );
	name = "scaff";
	index = 0;
	options[index++] = name;
	options[index++] = "-g";
	options[index++] = outfile;

	if ( getF )
	{
		options[index++] = "-F";
	}

	if ( getP )
	{
		options[index++] = "-p";
		sprintf ( ncpu_s, "%d", ncpu );
		options[index++] = ncpu_s;
	}

	if ( getL )
	{
		options[index++] = "-L";
		sprintf ( len_s, "%d", cutoff_len );
		options[index++] = len_s;
	}

	if ( getG )
	{
		options[index++] = "-G";
		sprintf ( gap_diff_s, "%d", gap_diff );
		options[index++] = gap_diff_s;
	}

	if ( getu )
	{
		options[index++] = "-u";
	}

	if ( getc )
	{
		options[index++] = "-c";
		sprintf ( min_cvg_s, "%f", min_cvg );
		options[index++] = min_cvg_s;
	}

	if ( getC )
	{
		options[index++] = "-C";
		sprintf ( max_cvg_s, "%f", max_cvg );
		options[index++] = max_cvg_s;
	}

	if ( getb )
	{
		options[index++] = "-b";
		sprintf ( insert_size_bound_s, "%f", insert_size_bound );
		options[index++] = insert_size_bound_s;
	}

	if ( getB )
	{
		options[index++] = "-B";
		sprintf ( bubble_coverage_s, "%f", bubble_coverage );
		options[index++] = bubble_coverage_s;
	}

	if ( getN )
	{
		options[index++] = "-N";
		sprintf ( genome_size_s, "%d", genome_size );
		options[index++] = genome_size_s;
	}

	if ( getw )
	{
		options[index++] = "-w";
	}

	if ( getV )
	{
		options[index++] = "-V";
	}

	/*
	for (i = 0; i < index; i++)
	{
	    fprintf (stderr,"%s ", options[i]);
	}

	fprintf (stderr,"\n");
	*/
	call_scaffold ( index, options );
	time ( &stop_t );
	fprintf ( stderr, "Time for the whole pipeline: %dm.\n", ( int ) ( stop_t - start_t ) / 60 );
}

static void display_all_usage ()
{
	//  fprintf (stderr,"\nSOAPdenovo all -s configFile -o outputGraph [-R -f -F -u -w] [-K kmer -p n_cpu -a initMemoryAssumption -d KmerFreqCutOff -D EdgeCovCutoff -M mergeLevel -k kmer_R2C, -G gapLenDiff -L minContigLen -c minContigCvg -C maxContigCvg -b insertSizeUpperBound -B bubbleCoverage -N genomeSize]\n");
	fprintf ( stderr, "\nSOAPdenovo all -s configFile -o outputGraph [-R -F -u -w] [-K kmer -p n_cpu -a initMemoryAssumption -d KmerFreqCutOff -D EdgeCovCutoff -M mergeLevel -k kmer_R2C, -G gapLenDiff -L minContigLen -c minContigCvg -C maxContigCvg -b insertSizeUpperBound -B bubbleCoverage -N genomeSize]\n" );
	fprintf ( stderr, "  -s <string>    configFile: the config file of solexa reads\n" );
	fprintf ( stderr, "  -o <string>    outputGraph: prefix of output graph file name\n" );
#ifdef MER127
	fprintf ( stderr, "  -K <int>       kmer(min 13, max 127): kmer size, [23]\n" );
#else
	fprintf ( stderr, "  -K <int>       kmer(min 13, max 63): kmer size, [23]\n" );
#endif
	fprintf ( stderr, "  -p <int>       n_cpu: number of cpu for use, [8]\n" );
	fprintf ( stderr, "  -a <int>       initMemoryAssumption: memory assumption initialized to avoid further reallocation, unit G, [0]\n" );
	fprintf ( stderr, "  -d <int>       kmerFreqCutoff: kmers with frequency no larger than KmerFreqCutoff will be deleted, [0]\n" );
	fprintf ( stderr, "  -R (optional)  resolve repeats by reads, [NO]\n" );
	fprintf ( stderr, "  -D <int>       edgeCovCutoff: edges with coverage no larger than EdgeCovCutoff will be deleted, [1]\n" );
	fprintf ( stderr, "  -M <int>       mergeLevel(min 0, max 3): the strength of merging similar sequences during contiging, [1]\n" );
	fprintf ( stderr, "  -e <int>       arcWeight: two edges, between which the arc's weight is larger than arcWeight, will be linerized, [0]\n" );
#ifdef MER127
	fprintf ( stderr, "  -m <int>       maxKmer (max 127): maximum kmer size used for multi-kmer, [NO]\n" );
#else
	fprintf ( stderr, "  -m <int>       maxKmer (max 63): maximum kmer size used for multi-kmer, [NO]\n" );
#endif
	fprintf ( stderr, "  -E (optional)  merge clean bubble before iterate, works only if -M is set when using multi-kmer, [NO]\n" );
	//  printf ("  -O (optional)\toutput contig of each kmer when iterating\n");
	//  fprintf (stderr,"  -f (optional)  output gap related reads in map step for using SRkgf to fill gaps, [NO]\n");
#ifdef MER127
	fprintf ( stderr, "  -k <int>       kmer_R2C(min 13, max 127): kmer size used for mapping reads to contigs, [K]\n" );
#else
	fprintf ( stderr, "  -k <int>       kmer_R2C(min 13, max 63): kmer size used for mapping reads to contigs, [K]\n" );
#endif
	fprintf ( stderr, "  -F (optional)  fill gaps in scaffolds, [NO]\n" );
	fprintf ( stderr, "  -u (optional)  un-mask contigs with high/low coverage before scaffolding, [mask]\n" );
	fprintf ( stderr, "  -w (optional)  keep contigs weakly connected to other contigs in scaffold, [NO]\n" );
	fprintf ( stderr, "  -G <int>       gapLenDiff: allowed length difference between estimated and filled gap, [50]\n" );
	fprintf ( stderr, "  -L <int>       minContigLen: shortest contig for scaffolding, [K+2]\n" );
	fprintf ( stderr, "  -c <float>     minContigCvg: minimum contig coverage (c*avgCvg), contigs shorter than 100bp with coverage smaller than c*avgCvg will be masked before scaffolding unless -u is set, [0.1]\n" );
	fprintf ( stderr, "  -C <float>     maxContigCvg: maximum contig coverage (C*avgCvg), contigs with coverage larger than C*avgCvg or contigs shorter than 100bp with coverage larger than 0.8*C*avgCvg will be masked before scaffolding unless -u is set, [2]\n" );
	fprintf ( stderr, "  -b <float>     insertSizeUpperBound: (b*avg_ins) will be used as upper bound of insert size for large insert size ( > 1000) when handling pair-end connections between contigs if b is set to larger than 1, [1.5]\n" );
	fprintf ( stderr, "  -B <float>     bubbleCoverage: remove contig with lower cvoerage in bubble structure if both contigs' coverage are smaller than bubbleCoverage*avgCvg, [0.6]\n" );
	fprintf ( stderr, "  -N <int>       genomeSize: genome size for statistics, [0]\n" );
	fprintf ( stderr, "  -V (optional)  output information for Hawkeye to visualize the assembly, [NO]\n" );
}

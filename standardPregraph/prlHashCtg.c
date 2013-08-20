/*
 * prlHashCtg.c
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

//debugging variables
static long long * kmerCounter;

//buffer related varibles for chop kmer
static unsigned int read_c;
static char ** rcSeq;
static char * seqBuffer;
static int * lenBuffer;
static unsigned int * indexArray;
static unsigned int * seqBreakers;  //record sum length to indicate start pos on seqBuffer
static int * ctgIdArray;

//static Kmer *firstKmers;

//buffer related varibles for splay tree
static unsigned int buffer_size = 100000000;
static unsigned int seq_buffer_size;
static unsigned int max_read_c;
static volatile unsigned int kmer_c;
static Kmer * kmerBuffer;
static ubyte8 * hashBanBuffer;
static boolean * smallerBuffer;

static void singleKmer ( int t, KmerSet * kset, unsigned int seq_index, unsigned int pos );
static void chopKmer4read ( int t, int threadID );

static void threadRoutine ( void * para )
{
	PARAMETER * prm;
	unsigned int i;
	unsigned char id;
	prm = ( PARAMETER * ) para;
	id = prm->threadID;

	while ( 1 )
	{
		if ( * ( prm->selfSignal ) == 1 )
		{
			unsigned int seq_index = 0;
			unsigned int pos = 0;

			for ( i = 0; i < kmer_c; i++ )
			{
				if ( seq_index < read_c && indexArray[seq_index + 1] == i )
				{
					seq_index++;    // which sequence this kmer belongs to
					pos = 0;
				}

				if ( ( unsigned char ) ( hashBanBuffer[i] % thrd_num ) != id )
				{
					pos++;
					continue;
				}

				kmerCounter[id + 1]++;
				singleKmer ( i, KmerSets[id], seq_index, pos++ );
			}

			* ( prm->selfSignal ) = 0;
		}
		else if ( * ( prm->selfSignal ) == 2 )
		{
			for ( i = 0; i < read_c; i++ )
			{
				if ( i % thrd_num != id )
				{
					continue;
				}

				chopKmer4read ( i, id + 1 );
			}

			* ( prm->selfSignal ) = 0;
		}
		else if ( * ( prm->selfSignal ) == 3 )
		{
			* ( prm->selfSignal ) = 0;
			break;
		}

		usleep ( 1 );
	}
}

/*************************************************
Function:
    singleKmer
Description:
    Stores a kmer into kmerset and updates related status.
Input:
    1. t:               kmerBuffer index
    2. kset:            KmerSet
    3. seq_index:       contig id array index
    4. pos:         the postion of the kmer on contig
Output:
    None.
Return:
    None.
*************************************************/
static void singleKmer ( int t, KmerSet * kset, unsigned int seq_index, unsigned int pos )
{
	boolean flag;
	kmer_t * node;
	flag = put_kmerset ( kset, kmerBuffer[t], 4, 4, &node );

	if ( !flag )
	{
		if ( smallerBuffer[t] )
		{
			node->twin = 0;
		}
		else
		{
			node->twin = 1;
		};

		node->l_links = ctgIdArray[seq_index];

		node->r_links = pos;
	}
	else
	{
		node->deleted = 1;
	}
}

static void creatThrds ( pthread_t * threads, PARAMETER * paras )
{
	unsigned char i;
	int temp;

	for ( i = 0; i < thrd_num; i++ )
	{
		if ( ( temp = pthread_create ( &threads[i], NULL, ( void * ) threadRoutine, & ( paras[i] ) ) ) != 0 )
		{
			fprintf ( stderr, "Create threads failed.\n" );
			exit ( 1 );
		}
	}

	fprintf ( stderr, "%d thread(s) initialized.\n", thrd_num );
}

static void thread_wait ( pthread_t * threads )
{
	int i;

	for ( i = 0; i < thrd_num; i++ )
		if ( threads[i] != 0 )
		{
			pthread_join ( threads[i], NULL );
		}
}

static void chopKmer4read ( int t, int threadID )
{
	char * src_seq = seqBuffer + seqBreakers[t];
	char * bal_seq = rcSeq[threadID];
	int len_seq = lenBuffer[t];
	int j, bal_j;
	ubyte8 hash_ban, bal_hash_ban;
	Kmer word, bal_word;
	int index;
#ifdef MER127
	word.high1 = word.low1 = word.high2 = word.low2 = 0;

	for ( index = 0; index < overlaplen; index++ )
	{
		word = KmerLeftBitMoveBy2 ( word );
		word.low2 |= src_seq[index];
	}

#else
	word.high = word.low = 0;

	for ( index = 0; index < overlaplen; index++ )
	{
		word = KmerLeftBitMoveBy2 ( word );
		word.low |= src_seq[index];
	}

#endif
	reverseComplementSeq ( src_seq, len_seq, bal_seq );
	// complementary node
	bal_word = reverseComplement ( word, overlaplen );
	bal_j = len_seq - 0 - overlaplen;
	index = indexArray[t];

	if ( KmerSmaller ( word, bal_word ) )
	{
		hash_ban = hash_kmer ( word );
		kmerBuffer[index] = word;
		hashBanBuffer[index] = hash_ban;
		smallerBuffer[index++] = 1;
	}
	else
	{
		bal_hash_ban = hash_kmer ( bal_word );
		kmerBuffer[index] = bal_word;
		hashBanBuffer[index] = bal_hash_ban;
		smallerBuffer[index++] = 0;
	}

	for ( j = 1; j <= len_seq - overlaplen; j++ )
	{
		word = nextKmer ( word, src_seq[j - 1 + overlaplen] );
		bal_j = len_seq - j - overlaplen;
		bal_word = prevKmer ( bal_word, bal_seq[bal_j] );

		if ( KmerSmaller ( word, bal_word ) )
		{
			hash_ban = hash_kmer ( word );
			kmerBuffer[index] = word;
			hashBanBuffer[index] = hash_ban;
			smallerBuffer[index++] = 1;
		}
		else
		{
			// complementary node
			bal_hash_ban = hash_kmer ( bal_word );
			kmerBuffer[index] = bal_word;
			hashBanBuffer[index] = bal_hash_ban;
			smallerBuffer[index++] = 0;
		}
	}
}

static void sendWorkSignal ( unsigned char SIG, unsigned char * thrdSignals )
{
	int t;

	for ( t = 0; t < thrd_num; t++ )
	{
		thrdSignals[t + 1] = SIG;
	}

	while ( 1 )
	{
		usleep ( 10 );

		for ( t = 0; t < thrd_num; t++ )
			if ( thrdSignals[t + 1] )
			{
				break;
			}

		if ( t == thrd_num )
		{
			break;
		}
	}
}


/*************************************************
Function:
    getID
Description:
    Gets contig id.
Input:
    1.name      contig name
Output:
    None.
Return:
    None.
*************************************************/
static int getID ( char * name )
{
	if ( name[0] >= '0' && name[0] <= '9' )
	{
		return atoi ( & ( name[0] ) );
	}
	else
	{
		return 0;
	}
}


/*************************************************
Function:
    prlContig2nodes
Description:
    Chops contigs not shorter than 'len_cut'  into kmers:
    1. if the kmer occured twice or more , mark it deleted;
    2. 'node->l_links' stores the contig id and 'node->r_links' stores the position the kmer on contig.
Input:
    1. grapfile:        the prefix of contig file
    2. len_cut:     the minimum length of contig
Output:
    None.
Return:
    True always.
*************************************************/
boolean prlContig2nodes ( char * grapfile, int len_cut )
{
	long long i, num_seq;
	char name[256], *next_name;
	FILE * fp;
	pthread_t threads[thrd_num];
	time_t start_t, stop_t;
	unsigned char thrdSignal[thrd_num + 1];
	PARAMETER paras[thrd_num];
	int maxCtgLen, minCtgLen, nameLen;
	unsigned int lenSum, contigId;
	//init
	WORDFILTER = createFilter ( overlaplen );
	time ( &start_t );
	sprintf ( name, "%s.contig", grapfile );
	fp = ckopen ( name, "r" );
	maxCtgLen = nameLen = 10;
	minCtgLen = 1000;
	num_seq = readseqpar ( &maxCtgLen, &minCtgLen, &nameLen, fp );
	fprintf ( stderr, "\n%lld contig(s), maximum sequence length %d, minimum sequence length %d, maximum name length %d.\n", num_seq, maxCtgLen, minCtgLen, nameLen );
	maxReadLen = maxCtgLen;
	fclose ( fp );
	time ( &stop_t );
	fprintf ( stderr, "Time spent on parsing contigs file: %ds.\n", ( int ) ( stop_t - start_t ) );
	next_name = ( char * ) ckalloc ( ( maxNameLen + 1 ) * sizeof ( char ) );
	// extract all the EDONs
	seq_buffer_size = buffer_size * 2;
	max_read_c = seq_buffer_size / 20;
	kmerBuffer = ( Kmer * ) ckalloc ( buffer_size * sizeof ( Kmer ) );
	hashBanBuffer = ( ubyte8 * ) ckalloc ( buffer_size * sizeof ( ubyte8 ) );
	smallerBuffer = ( boolean * ) ckalloc ( buffer_size * sizeof ( boolean ) );
	seqBuffer = ( char * ) ckalloc ( seq_buffer_size * sizeof ( char ) );
	lenBuffer = ( int * ) ckalloc ( max_read_c * sizeof ( int ) );
	indexArray = ( unsigned int * ) ckalloc ( ( max_read_c + 1 ) * sizeof ( unsigned int ) );
	seqBreakers = ( unsigned int * ) ckalloc ( ( max_read_c + 1 ) * sizeof ( unsigned int ) );
	ctgIdArray = ( int * ) ckalloc ( max_read_c * sizeof ( int ) );
	fp = ckopen ( name, "r" );
	rcSeq = ( char ** ) ckalloc ( ( thrd_num + 1 ) * sizeof ( char * ) );

	if ( 1 )
	{
		kmerCounter = ( long long * ) ckalloc ( ( thrd_num + 1 ) * sizeof ( long long ) );
		KmerSets = ( KmerSet ** ) ckalloc ( thrd_num * sizeof ( KmerSet * ) );

		for ( i = 0; i < thrd_num; i++ )
		{
			KmerSets[i] = init_kmerset ( 1024, 0.77f );
			thrdSignal[i + 1] = 0;
			paras[i].threadID = i;
			paras[i].mainSignal = &thrdSignal[0];
			paras[i].selfSignal = &thrdSignal[i + 1];
			kmerCounter[i + 1] = 0;
			rcSeq[i + 1] = ( char * ) ckalloc ( maxCtgLen * sizeof ( char ) );
		}

		creatThrds ( threads, paras );
	}

	kmer_c = thrdSignal[0] = kmerCounter[0] = 0;
	time ( &start_t );
	read_c = lenSum = i = seqBreakers[0] = indexArray[0] = 0;
	readseq1by1 ( seqBuffer + seqBreakers[read_c], next_name, & ( lenBuffer[read_c] ), fp, -1 );

	while ( !feof ( fp ) )
	{
		contigId = getID ( next_name );
		readseq1by1 ( seqBuffer + seqBreakers[read_c], next_name, & ( lenBuffer[read_c] ), fp, 1 );

		if ( ( ++i ) % 10000000 == 0 )
		{
			fprintf ( stderr, "--- %lldth contig(s).\n", i );
		}

		if ( lenBuffer[read_c] < overlaplen + 1 || lenBuffer[read_c] < len_cut )
		{
			contigId = getID ( next_name );
			continue;
		}

		ctgIdArray[read_c] = contigId > 0 ? contigId : i;
		lenSum += lenBuffer[read_c];
		kmer_c += lenBuffer[read_c] - overlaplen + 1;
		read_c++;
		seqBreakers[read_c] = lenSum;
		indexArray[read_c] = kmer_c;

		if ( read_c == max_read_c || ( lenSum + maxCtgLen ) > seq_buffer_size || ( kmer_c + maxCtgLen - overlaplen + 1 ) > buffer_size )
		{
			kmerCounter[0] += kmer_c;
			sendWorkSignal ( 2, thrdSignal ); //chopKmer4read
			sendWorkSignal ( 1, thrdSignal ); //singleKmer
			kmer_c = read_c = lenSum = 0;
		}
	}

	if ( read_c )
	{
		kmerCounter[0] += kmer_c;
		sendWorkSignal ( 2, thrdSignal ); //chopKmer4read
		sendWorkSignal ( 1, thrdSignal ); //singleKmer
	}

	sendWorkSignal ( 3, thrdSignal ); //stop threads
	thread_wait ( threads );
	time ( &stop_t );
	fprintf ( stderr, "Time spent on hashing contigs: %ds.\n", ( int ) ( stop_t - start_t ) );

	if ( 1 )
	{
		unsigned long long alloCounter = 0;
		unsigned long long allKmerCounter = 0;

		for ( i = 0; i < thrd_num; i++ )
		{
			alloCounter += count_kmerset ( ( KmerSets[i] ) );
			allKmerCounter += kmerCounter[i + 1];
			free ( ( void * ) rcSeq[i + 1] );
		}

		fprintf ( stderr, "%lli node(s) allocated, %lli kmer(s) in contigs, %lli kmer(s) processed.\n", alloCounter, kmerCounter[0], allKmerCounter );
	}

	free ( ( void * ) rcSeq );
	free ( ( void * ) kmerCounter );
	free ( ( void * ) seqBuffer );
	free ( ( void * ) lenBuffer );
	free ( ( void * ) indexArray );
	free ( ( void * ) seqBreakers );
	free ( ( void * ) ctgIdArray );
	free ( ( void * ) kmerBuffer );
	free ( ( void * ) hashBanBuffer );
	free ( ( void * ) smallerBuffer );
	free ( ( void * ) next_name );
	fclose ( fp );
	return 1;
}

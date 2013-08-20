/*
 * prlHashReads.c
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
static long long * tips;
static long long * kmerCounter; //kmer number for each KmerSet[thread id]

static long long ** kmerFreq;

//buffer related varibles for chop kmer
static int read_c;
static char ** rcSeq;       //sequence pointer for each thread
static char ** seqBuffer;   //read buffer
static int * lenBuffer;     //read length buffer
static int * indexArray; // the read's begin kmer 's kmer_c

//buffer related varibles for splay tree
static int buffer_size = 100000000;     //buffer size for kmerBuffer...
static volatile int kmer_c;             // kmer count number
static Kmer * kmerBuffer;               //kmer buffer array
static ubyte8 * hashBanBuffer;          //the buffered hash value for 'kmerBuffer'
static char * nextcBuffer, *prevcBuffer; //next char buffer , previous char buffer for 'kmerBuffer'

static struct aiocb aio1;
static struct aiocb aio2;
static char * aioBuffer1;
static char * aioBuffer2;
static char * readBuffer1;
static char * readBuffer2;

static void thread_mark ( KmerSet * set, unsigned char thrdID );
static void Mark1in1outNode ( unsigned char * thrdSignal );
static void thread_delow ( KmerSet * set, unsigned char thrdID );
static void deLowCov ( unsigned char * thrdSignal );

static void singleKmer ( int t, KmerSet * kset );
static void chopKmer4read ( int t, int threadID );

static void freqStat ( char * outfile );

static void threadRoutine ( void * para )
{
	PARAMETER * prm;
	int i;
	unsigned char id;
	prm = ( PARAMETER * ) para;
	id = prm->threadID;

	//printf("%dth thread with threadID %d, hash_table %p\n",id,prm.threadID,prm.hash_table);
	while ( 1 )
	{
		if ( * ( prm->selfSignal ) == 1 )
		{
			for ( i = 0; i < kmer_c; i++ )
			{
				//if((unsigned char)(magic_seq(hashBanBuffer[i])%thrd_num)!=id)
				//if((kmerBuffer[i]%thrd_num)!=id)
				if ( ( hashBanBuffer[i] % thrd_num ) != id )
				{
					continue;
				}

				kmerCounter[id + 1]++;
				singleKmer ( i, KmerSets[id] );
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
		else if ( * ( prm->selfSignal ) == 4 )
		{
			thread_mark ( KmerSets[id], id );
			* ( prm->selfSignal ) = 0;
		}
		else if ( * ( prm->selfSignal ) == 5 )
		{
			thread_delow ( KmerSets[id], id );
			* ( prm->selfSignal ) = 0;
		}

		usleep ( 1 );
	}
}

static void singleKmer ( int t, KmerSet * kset )
{
	kmer_t * pos;
	put_kmerset ( kset, kmerBuffer[t], prevcBuffer[t], nextcBuffer[t], &pos );
}

static void creatThrds ( pthread_t * threads, PARAMETER * paras )
{
	unsigned char i;
	int temp;

	for ( i = 0; i < thrd_num; i++ )
	{
		//printf("to create %dth thread\n",(*(char *)&(threadID[i])));
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
	char * src_seq = seqBuffer[t];
	char * bal_seq = rcSeq[threadID];
	int len_seq = lenBuffer[t];
	int j, bal_j;
	ubyte8 hash_ban, bal_hash_ban;
	Kmer word, bal_word;
	int index;
	char InvalidCh = 4;
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
	bal_j = len_seq  - overlaplen;
	index = indexArray[t];

	if ( KmerSmaller ( word, bal_word ) )
	{
		hash_ban = hash_kmer ( word );
		hashBanBuffer[index] = hash_ban;
		kmerBuffer[index] = word;
		prevcBuffer[index] = InvalidCh;
		nextcBuffer[index++] = src_seq[0 + overlaplen];
	}
	else
	{
		bal_hash_ban = hash_kmer ( bal_word );
		hashBanBuffer[index] = bal_hash_ban;
		kmerBuffer[index] = bal_word;
		prevcBuffer[index] = bal_seq[bal_j - 1];
		nextcBuffer[index++] = InvalidCh;
	}

	for ( j = 1; j <= len_seq - overlaplen; j++ )
	{
		word = nextKmer ( word, src_seq[j - 1 + overlaplen] );
		bal_j = len_seq - j - overlaplen;
		bal_word = prevKmer ( bal_word, bal_seq[bal_j] );

		if ( KmerSmaller ( word, bal_word ) )
		{
			hash_ban = hash_kmer ( word );
			hashBanBuffer[index] = hash_ban;
			kmerBuffer[index] = word;
			prevcBuffer[index] = src_seq[j - 1];

			if ( j < len_seq - overlaplen )
			{
				nextcBuffer[index++] = src_seq[j + overlaplen];
			}
			else
			{
				nextcBuffer[index++] = InvalidCh;
			}

			//printf("%dth: %p with %p\n",kmer_c-1,word,hashBanBuffer[kmer_c-1]);
		}
		else
		{
			// complementary node
			bal_hash_ban = hash_kmer ( bal_word );
			hashBanBuffer[index] = bal_hash_ban;
			kmerBuffer[index] = bal_word;

			if ( bal_j > 0 )
			{
				prevcBuffer[index] = bal_seq[bal_j - 1];
			}
			else
			{
				prevcBuffer[index] = InvalidCh;
			}

			nextcBuffer[index++] = bal_seq[bal_j + overlaplen];
			//printf("%dth: %p with %p\n",kmer_c-1,bal_word,hashBanBuffer[kmer_c-1]);
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
    prlRead2HashTable
Description:
    1. Imports the reads from the lib file one by one.
    2. Chops the reads into kmers and store them in KmerSets.
    3. Removes the kmers with low coverage.
    4. Marks the linear kmers.
    5. Counts the kmer frequences.
Input:
    1. libfile :            the reads config file
    2. outfile :        the output file prefix
Output:
    None.
Return:
    1 if exits normally.
*************************************************/
boolean prlRead2HashTable ( char * libfile, char * outfile )
{
	char * cach1;
	char * cach2;
	unsigned char asm_ctg = 1;
	long long i;
	char * next_name, name[256];
	FILE * fo;
	time_t start_t, stop_t;
	int maxReadNum;
	int libNo;
	pthread_t threads[thrd_num];
	unsigned char thrdSignal[thrd_num + 1];
	PARAMETER paras[thrd_num];
	boolean flag, pairs = 0;
	WORDFILTER = createFilter ( overlaplen );
	maxReadLen = 0;
	maxNameLen = 256;
	scan_libInfo ( libfile );
	alloc_pe_mem ( num_libs );

	if ( !maxReadLen )
	{
		maxReadLen = 100;
	}

	if ( gLineLen < maxReadLen )
	{
		gStr = ( char * ) ckalloc ( ( maxReadLen + 1 ) * sizeof ( char ) );
	}

	//init
	maxReadLen4all = maxReadLen;
	fprintf ( stderr, "In %s, %d lib(s), maximum read length %d, maximum name length %d.\n\n", libfile, num_libs, maxReadLen, maxNameLen );
	next_name = ( char * ) ckalloc ( ( maxNameLen + 1 ) * sizeof ( char ) );
	kmerBuffer = ( Kmer * ) ckalloc ( buffer_size * sizeof ( Kmer ) );
	hashBanBuffer = ( ubyte8 * ) ckalloc ( buffer_size * sizeof ( ubyte8 ) );
	prevcBuffer = ( char * ) ckalloc ( buffer_size * sizeof ( char ) );
	nextcBuffer = ( char * ) ckalloc ( buffer_size * sizeof ( char ) );
	maxReadNum = buffer_size / ( maxReadLen - overlaplen + 1 );
	//printf("buffer size %d, max read len %d, max read num %d\n",buffer_size,maxReadLen,maxReadNum);
	int maxAIOSize = 32768;
	aioBuffer1 = ( char * ) ckalloc ( ( maxAIOSize ) * sizeof ( char ) );
	aioBuffer2 = ( char * ) ckalloc ( ( maxAIOSize ) * sizeof ( char ) );
	readBuffer1 = ( char * ) ckalloc ( ( maxAIOSize + ( maxReadLen * 4 + 1024 ) ) * sizeof ( char ) ); //(char *)ckalloc(maxAIOSize*sizeof(char)); //1024
	readBuffer2 = ( char * ) ckalloc ( ( maxAIOSize + ( maxReadLen * 4 + 1024 ) ) * sizeof ( char ) ); //1024
	cach1 = ( char * ) ckalloc ( ( maxReadLen * 4 + 1024 ) * sizeof ( char ) ); //1024
	cach2 = ( char * ) ckalloc ( ( maxReadLen * 4 + 1024 ) * sizeof ( char ) ); //1024
	memset ( cach1, '\0', ( maxReadLen * 4 + 1024 ) ); //1024
	memset ( cach2, '\0', ( maxReadLen * 4 + 1024 ) ); //1024
	seqBuffer = ( char ** ) ckalloc ( maxReadNum * sizeof ( char * ) );
	lenBuffer = ( int * ) ckalloc ( maxReadNum * sizeof ( int ) );
	indexArray = ( int * ) ckalloc ( maxReadNum * sizeof ( int ) );

	for ( i = 0; i < maxReadNum; i++ )
	{
		seqBuffer[i] = ( char * ) ckalloc ( maxReadLen * sizeof ( char ) );
	}

	rcSeq = ( char ** ) ckalloc ( ( thrd_num + 1 ) * sizeof ( char * ) );

	if ( 1 )
	{
		kmerCounter = ( long long * ) ckalloc ( ( thrd_num + 1 ) * sizeof ( long long ) );
		KmerSets = ( KmerSet ** ) ckalloc ( thrd_num * sizeof ( KmerSet * ) );
		ubyte8 init_size = 1024;
		ubyte8 k = 0;

		if ( initKmerSetSize )
		{
#ifdef MER127
			init_size = ( ubyte8 ) ( ( double ) initKmerSetSize * 1024.0f * 1024.0f * 1024.0f / ( double ) thrd_num / 40 );
#else
			init_size = ( ubyte8 ) ( ( double ) initKmerSetSize * 1024.0f * 1024.0f * 1024.0f / ( double ) thrd_num / 24 ); //is it true?
#endif

			do
			{
				++k;
			}
			while ( k * 0xFFFFFFLLU < init_size );
		}

		for ( i = 0; i < thrd_num; i++ )
		{
			//KmerSets[i] = init_kmerset(1024,0.77f);
			KmerSets[i] = init_kmerset ( ( ( initKmerSetSize ) ? ( k * 0xFFFFFFLLU ) : ( init_size ) ), 0.77f );
			thrdSignal[i + 1] = 0;
			paras[i].threadID = i;
			paras[i].mainSignal = &thrdSignal[0];
			paras[i].selfSignal = &thrdSignal[i + 1];
			kmerCounter[i + 1] = 0;
			rcSeq[i + 1] = ( char * ) ckalloc ( maxReadLen * sizeof ( char ) );
		}

		creatThrds ( threads, paras );
	}

	thrdSignal[0] = kmerCounter[0] = 0;
	time ( &start_t );
	kmer_c = n_solexa = read_c = i = libNo = readNumBack = gradsCounter = 0;

	while ( openNextFile ( &libNo, pairs, asm_ctg ) )
	{
		//read bam file
		if ( lib_array[libNo].curr_type == 4 )
		{
			int type = 0;   //deside the PE reads is good or bad

			while ( ( flag = read1seqInLibBam ( seqBuffer[read_c], next_name, & ( lenBuffer[read_c] ), &libNo, pairs, 1, &type ) ) != 0 )
			{
				if ( type == -1 ) //if the reads is bad, go back.
				{
					i--;

					if ( lenBuffer[read_c - 1] >= overlaplen + 1 )
					{
						kmer_c -= lenBuffer[read_c - 1] - overlaplen + 1;
						read_c--;
					}

					n_solexa -= 2;
					continue;
				}

				if ( ( ++i ) % 100000000 == 0 )
					{ fprintf ( stderr, "--- %lldth reads.\n", i ); }

				if ( lenBuffer[read_c] < 0 )
					{ fprintf ( stderr, "Read len %d.\n", lenBuffer[read_c] ); }

				if ( lenBuffer[read_c] < overlaplen + 1 )
					{ continue; }

				/*
				   if(lenBuffer[read_c]>70)
				   lenBuffer[read_c] = 50;
				   else if(lenBuffer[read_c]>40)
				   lenBuffer[read_c] = 40;
				 */
				indexArray[read_c] = kmer_c;
				kmer_c += lenBuffer[read_c] - overlaplen + 1;
				read_c++;

				if ( read_c == maxReadNum )
				{
					kmerCounter[0] += kmer_c;
					sendWorkSignal ( 2, thrdSignal ); //chopKmer4read
					sendWorkSignal ( 1, thrdSignal ); //singleKmer
					kmer_c = read_c = 0;
				}
			}
		}
		//read PE fasta or fastq
		else if ( lib_array[libNo].curr_type == 1 || lib_array[libNo].curr_type == 2 )
		{
			initAIO ( &aio1, aioBuffer1, fileno ( lib_array[libNo].fp1 ), maxAIOSize );
			initAIO ( &aio2, aioBuffer2, fileno ( lib_array[libNo].fp2 ), maxAIOSize );
			int offset1, offset2, flag1, flag2, rt1, rt2;
			offset1 = offset2 = 0;
			rt1 = aio_read ( &aio1 );
			rt2 = aio_read ( &aio2 );
			flag1 = AIORead ( &aio1, &offset1, readBuffer1, cach1, &rt1, lib_array[libNo].curr_type );
			flag2 = AIORead ( &aio2, &offset2, readBuffer2, cach2, &rt2, lib_array[libNo].curr_type );

			if ( flag1 && flag2 )
			{
				int start1, start2, turn;
				start1 = start2 = 0;
				turn = 1;

				while ( start1 < offset1 || start2 < offset2 )
				{
					if ( turn == 1 )
					{
						turn = 2;
						readseqInLib ( seqBuffer[read_c], next_name, & ( lenBuffer[read_c] ), readBuffer1, &start1, offset1, libNo );

						if ( ( ++i ) % 100000000 == 0 )
							{ fprintf ( stderr, "--- %lldth reads.\n", i ); }

						if ( lenBuffer[read_c] < 0 )
							{ fprintf ( stderr, "Read len %d.\n", lenBuffer[read_c] ); }

						if ( lenBuffer[read_c] < overlaplen + 1 )
						{
							if ( start1 >= offset1 )
							{
								start1 = 0;
								offset1 = 0;
								flag1 = AIORead ( &aio1, &offset1, readBuffer1, cach1, &rt1, lib_array[libNo].curr_type );
							}

							continue;
						}

						indexArray[read_c] = kmer_c;
						kmer_c += lenBuffer[read_c] - overlaplen + 1;
						read_c++;

						if ( start1 >= offset1 )
						{
							start1 = 0;
							offset1 = 0;
							flag1 = AIORead ( &aio1, &offset1, readBuffer1, cach1, &rt1, lib_array[libNo].curr_type );
						}

						if ( read_c == maxReadNum )
						{
							kmerCounter[0] += kmer_c;
							sendWorkSignal ( 2, thrdSignal );   //chopKmer4read
							sendWorkSignal ( 1, thrdSignal );   //singleKmer
							kmer_c = read_c = 0;
						}

						continue;
					}

					if ( turn == 2 )
					{
						turn = 1;
						readseqInLib ( seqBuffer[read_c], next_name, & ( lenBuffer[read_c] ), readBuffer2, &start2, offset2, libNo );

						if ( ( ++i ) % 100000000 == 0 )
							{ fprintf ( stderr, "--- %lldth reads.\n", i ); }

						if ( lenBuffer[read_c] < 0 )
							{ fprintf ( stderr, "Read len %d.\n", lenBuffer[read_c] ); }

						if ( lenBuffer[read_c] < overlaplen + 1 )
						{
							if ( ( flag2 == 2 ) && ( start2 >= offset2 ) )
								{ break; }

							if ( start2 >= offset2 )
							{
								start2 = 0;
								offset2 = 0;
								flag2 = AIORead ( &aio2, &offset2, readBuffer2, cach2, &rt2, lib_array[libNo].curr_type );
							}

							continue;
						}

						indexArray[read_c] = kmer_c;
						kmer_c += lenBuffer[read_c] - overlaplen + 1;
						read_c++;

						if ( ( flag2 == 2 ) && ( start2 >= offset2 ) )
							{ break; }

						if ( start2 >= offset2 )
						{
							start2 = 0;
							offset2 = 0;
							flag2 = AIORead ( &aio2, &offset2, readBuffer2, cach2, &rt2, lib_array[libNo].curr_type );
						}

						if ( read_c == maxReadNum )
						{
							kmerCounter[0] += kmer_c;
							sendWorkSignal ( 2, thrdSignal );   //chopKmer4read
							sendWorkSignal ( 1, thrdSignal );   //singleKmer
							kmer_c = read_c = 0;
						}

						continue;
					}
				}
			}
			else
			{
				fprintf(stderr, "Error: aio_read error.\n");
			}
		}
		//read single fasta, single fastq and PE fasta in one file
		else
		{
			initAIO ( &aio1, aioBuffer1, fileno ( lib_array[libNo].fp1 ), maxAIOSize );
			int offset, flag1, rt;
			offset = 0;
			rt = aio_read ( &aio1 );

			while ( ( flag1 = AIORead ( &aio1, &offset, readBuffer1, cach1, &rt, lib_array[libNo].curr_type ) ) )
			{
				int start = 0;

				while ( start < offset )
				{
					readseqInLib ( seqBuffer[read_c], next_name, & ( lenBuffer[read_c] ), readBuffer1, &start, offset, libNo );

					if ( ( ++i ) % 100000000 == 0 )
						{ fprintf ( stderr, "--- %lldth reads.\n", i ); }

					if ( lenBuffer[read_c] < 0 )
						{ fprintf ( stderr, "Read len %d.\n", lenBuffer[read_c] ); }

					if ( lenBuffer[read_c] < overlaplen + 1 )
						{ continue; }

					indexArray[read_c] = kmer_c;
					kmer_c += lenBuffer[read_c] - overlaplen + 1;
					read_c++;
				}

				if ( read_c > maxReadNum - 1024 )
				{
					kmerCounter[0] += kmer_c;
					sendWorkSignal ( 2, thrdSignal );   //chopKmer4read
					sendWorkSignal ( 1, thrdSignal );   //singleKmer
					kmer_c = read_c = 0;
				}

				if ( flag1 == 2 )
					{ break; }
			}
		}
	}

	if ( read_c )
	{
		kmerCounter[0] += kmer_c;
		sendWorkSignal ( 2, thrdSignal );   //chopKmer4read
		sendWorkSignal ( 1, thrdSignal );   //singleKmer
	}

	time ( &stop_t );
	fprintf ( stderr, "Time spent on hashing reads: %ds, %lld read(s) processed.\n", ( int ) ( stop_t - start_t ), i );

	//record insert size info
	if ( pairs )
	{
		if ( gradsCounter )
			{ fprintf ( stderr, "%d pe insert size, the largest boundary is %lld.\n\n", gradsCounter, pes[gradsCounter - 1].PE_bound ); }
		else
		{
			fprintf ( stderr, "No paired reads found.\n" );
		}

		sprintf ( name, "%s.peGrads", outfile );
		fo = ckopen ( name, "w" );
		fprintf ( fo, "grads&num: %d\t%lld\n", gradsCounter, n_solexa );

		for ( i = 0; i < gradsCounter; i++ )
		{
			fprintf ( fo, "%d\t%lld\t%d\n", pes[i].insertS, pes[i].PE_bound, pes[i].rank );
		}

		fclose ( fo );
	}

	free_pe_mem ();
	free_libs ();

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

		fprintf ( stderr, "%lli node(s) allocated, %lli kmer(s) in reads, %lli kmer(s) processed.\n", alloCounter, kmerCounter[0], allKmerCounter );
	}

	free ( ( void * ) rcSeq );
	free ( ( void * ) kmerCounter );

	for ( i = 0; i < maxReadNum; i++ )
	{
		free ( ( void * ) seqBuffer[i] );
	}

	free ( ( void * ) seqBuffer );
	free ( ( void * ) lenBuffer );
	free ( ( void * ) indexArray );
	free ( ( void * ) kmerBuffer );
	free ( ( void * ) hashBanBuffer );
	free ( ( void * ) nextcBuffer );
	free ( ( void * ) prevcBuffer );
	free ( ( void * ) next_name );
	free ( ( void * ) aioBuffer1 );
	free ( ( void * ) aioBuffer2 );
	free ( ( void * ) readBuffer1 );
	free ( ( void * ) readBuffer2 );
	free ( ( void * ) cach1 );
	free ( ( void * ) cach2 );
	fprintf ( stderr, "done hashing nodes\n" );

	if ( deLowKmer )
	{
		time ( &start_t );
		deLowCov ( thrdSignal );
		time ( &stop_t );
		fprintf ( stderr, "Time spent on delowcvgNode: %ds.\n", ( int ) ( stop_t - start_t ) );
	}

	time ( &start_t );
	Mark1in1outNode ( thrdSignal );
	freqStat ( outfile );
	time ( &stop_t );
	fprintf ( stderr, "Time spent on marking linear nodes: %ds.\n", ( int ) ( stop_t - start_t ) );
	sendWorkSignal ( 3, thrdSignal );   //exit
	thread_wait ( threads );
	return 1;
}

void initAIO ( struct aiocb * aio, char * buf, int fd, int size )
{
	bzero ( aio, sizeof ( struct aiocb ) );
	aio->aio_buf = ( void * ) buf;
	aio->aio_fildes = fd;
	aio->aio_nbytes = size;
	aio->aio_offset = 0;
}

int AIORead ( struct aiocb * mycb, int * offset, char * buf, char * cach, int * rt, int curr_type )
{
	int i, i2, i3, j;
	int num;
	size_t mode, get, max_list;

	//      rt = aio_read(mycb);
	if ( *rt == 0 )
	{
		struct aiocb * aiocb_list[1];
		aiocb_list[0] = mycb;
		max_list = 1;

		while ( 1 )
		{
			mode = aio_suspend ( ( const struct aiocb * const * ) aiocb_list, max_list, NULL );

			if ( mode == -1 )
			{
				if ( errno != EAGAIN && errno != EINTR )
				{
					fprintf ( stderr, "Error:%s.\n", errno );
					return 0;
				}
				else
					{ continue; }
			}
			else
			{
				//while(aio_error(mycb) == EINPROGRESS);
				get = aio_return ( mycb );
				j = strlen ( cach );

				if ( get > 0 )
				{
					char * temp = ( char * ) ( ( *mycb ).aio_buf );

					if ( ( get % 32768 ) != 0 )
					{
						strcpy ( buf, cach );
						memcpy ( &buf[j], temp, get );
						memset ( cach, '\0', j );
						//printf("%s",buf);
						*offset = j + get;
						return 2;
					}

					if ( ( curr_type == 2 ) || ( curr_type == 6 ) )
					{
						num = 0;

						for ( i = get - 1; ( temp[i] != '@' ) || ( temp[i - 1] != '\n' ); i-- )
						{
							if ( temp[i] == '\n' ) {num++;}
						}

						if ( num <= 1 )
						{
							for ( i2 = i - 2; temp[i2] != '\n'; i2-- ) { ; }

							if ( temp[i2 + 1] == '+' )
							{
								for ( i2 = i2 - 1; temp[i2] != '\n'; i2-- ) { ; }

								if ( temp[i2 + 1] != '+' ) {for ( i = i2 - 1; ( temp[i] != '@' ) || ( temp[i - 1] != '\n' ); i-- ) { ; }}
							}
						}
					}
					else if ( ( curr_type == 1 ) || ( curr_type == 3 ) || ( curr_type == 5 ) )
						for ( i = get - 1; temp[i] != '>'; i-- ) { ; }

					//for (i = get - 1; temp[i] != '>' && temp[i] != '@'; i--) ;
					strcpy ( buf, cach );
					memcpy ( &buf[j], temp, i );
					//printf("%s",buf);
					*offset = i + j;
					memset ( cach, '\0', j );
					memcpy ( cach, &temp[i], get - i );
					( *mycb ).aio_offset += get;
					*rt = aio_read ( mycb );
					return 1;
				}
				else
				{
					fprintf(stderr, "Error: aio_return error.\n");
				}
				/*else
				   {
				   char *temp = (char *)((*mycb).aio_buf);
				   strcpy(buf,cach);
				   strcpy(&buf[j],temp);
				   *offset = j + get;
				   return 2;
				   } */
			}
		}
	}
	else
	{
		fprintf(stderr, "Error: (*rt != 0) in AIORead.\n");
	}

	return 0;
}

boolean openNextFile ( int * libNo, boolean pairs, unsigned char asm_ctg )
{
	int i = *libNo;
	int prevLib = i;

	if ( lib_array[i].fp1 )
		{ closeFp1InLab ( i ); }

	if ( lib_array[i].fp2 )
		{ closeFp2InLab ( i ); }

	*libNo = nextValidIndex ( i, pairs, asm_ctg );
	i = *libNo;

	if ( lib_array[i].rd_len_cutoff > 0 )
		{ maxReadLen = lib_array[i].rd_len_cutoff < maxReadLen4all ? lib_array[i].rd_len_cutoff : maxReadLen4all; }
	else
		{ maxReadLen = maxReadLen4all; }

	//record insert size info
	//printf("from lib %d to %d, read %lld to %ld\n",prevLib,i,readNumBack,n_solexa);
	if ( pairs && i != prevLib )
	{
		if ( readNumBack < n_solexa )
		{
			pes[gradsCounter].PE_bound = n_solexa;
			pes[gradsCounter].rank = lib_array[prevLib].rank;
			pes[gradsCounter].pair_num_cut = lib_array[prevLib].pair_num_cut;
			pes[gradsCounter++].insertS = lib_array[prevLib].avg_ins;
			readNumBack = n_solexa;
		}
	}

	if ( i >= num_libs )
		{ return 0; }

	openFileInLib ( i );
	return 1;
}

static void thread_delow ( KmerSet * set, unsigned char thrdID )
{
	int i, in_num, out_num, cvgSingle;
	int l_cvg, r_cvg;
	kmer_t * rs;
	set->iter_ptr = 0;

	while ( set->iter_ptr < set->size )
	{
		if ( !is_kmer_entity_null ( set->flags, set->iter_ptr ) )
		{
			in_num = out_num = l_cvg = r_cvg = 0;
			rs = set->array + set->iter_ptr;

			for ( i = 0; i < 4; i++ )
			{
				cvgSingle = get_kmer_left_cov ( *rs, i );

				if ( cvgSingle > 0 && cvgSingle <= deLowKmer )
				{
					set_kmer_left_cov ( *rs, i, 0 );
				}

				cvgSingle = get_kmer_right_cov ( *rs, i );

				if ( cvgSingle > 0 && cvgSingle <= deLowKmer )
				{
					set_kmer_right_cov ( *rs, i, 0 );
				}
			}

			if ( rs->l_links == 0 && rs->r_links == 0 )
			{
				rs->deleted = 1;
				tips[thrdID]++;
			}
		}

		set->iter_ptr++;
	}

	//printf("%lld single nodes, %lld linear\n",counter,tips[thrdID]);
}

static void deLowCov ( unsigned char * thrdSignal )
{
	int i;
	long long counter = 0;
	tips = ( long long * ) ckalloc ( thrd_num * sizeof ( long long ) );

	for ( i = 0; i < thrd_num; i++ )
	{
		tips[i] = 0;
	}

	sendWorkSignal ( 5, thrdSignal ); //mark linear nodes

	for ( i = 0; i < thrd_num; i++ )
	{
		counter += tips[i];
	}

	free ( ( void * ) tips );
	fprintf ( stderr, "%lld kmer(s) removed.\n", counter );
}

static void thread_mark ( KmerSet * set, unsigned char thrdID )
{
	int i, in_num, out_num, cvgSingle;
	int l_cvg, r_cvg;
	kmer_t * rs;
	long long counter = 0;
	set->iter_ptr = 0;

	while ( set->iter_ptr < set->size )
	{
		if ( !is_kmer_entity_null ( set->flags, set->iter_ptr ) )
		{
			in_num = out_num = l_cvg = r_cvg = 0;
			rs = set->array + set->iter_ptr;

			for ( i = 0; i < 4; i++ )
			{
				cvgSingle = get_kmer_left_cov ( *rs, i );

				if ( cvgSingle > 0 )
				{
					in_num++;
					l_cvg += cvgSingle;
				}

				cvgSingle = get_kmer_right_cov ( *rs, i );

				if ( cvgSingle > 0 )
				{
					out_num++;
					r_cvg += cvgSingle;
				}
			}

			if ( rs->single )
			{
				kmerFreq[thrdID][1]++;
				counter++;
			}
			else
			{
				kmerFreq[thrdID][ ( l_cvg > r_cvg ? l_cvg : r_cvg )]++;
			}

			if ( in_num == 1 && out_num == 1 )
			{
				rs->linear = 1;
				tips[thrdID]++;
			}
		}

		set->iter_ptr++;
	}

	//printf("%lld single nodes, %lld linear\n",counter,tips[thrdID]);
}

static void Mark1in1outNode ( unsigned char * thrdSignal )
{
	int i;
	long long counter = 0;
	tips = ( long long * ) ckalloc ( thrd_num * sizeof ( long long ) );
	kmerFreq = ( long long ** ) ckalloc ( thrd_num * sizeof ( long long * ) );

	for ( i = 0; i < thrd_num; i++ )
	{
		kmerFreq[i] = ( long long * ) ckalloc ( 257 * sizeof ( long long ) );
		memset ( kmerFreq[i], 0, 257 * sizeof ( long long ) );
		tips[i] = 0;
	}

	sendWorkSignal ( 4, thrdSignal ); //mark linear nodes

	for ( i = 0; i < thrd_num; i++ )
	{
		counter += tips[i];
	}

	free ( ( void * ) tips );
	fprintf ( stderr, "%lld linear node(s) marked.\n", counter );
}

static void freqStat ( char * outfile )
{
	FILE * fo;
	char name[256];
	int i, j;
	long long sum;
	sprintf ( name, "%s.kmerFreq", outfile );
	fo = ckopen ( name, "w" );

	for ( i = 1; i < 256; i++ )
	{
		sum = 0;

		for ( j = 0; j < thrd_num; j++ )
		{
			sum += kmerFreq[j][i];
		}

		fprintf ( fo, "%lld\n", sum );
	}

	for ( i = 0; i < thrd_num; i++ )
	{
		free ( ( void * ) kmerFreq[i] );
	}

	free ( ( void * ) kmerFreq );
	fclose ( fo );
}

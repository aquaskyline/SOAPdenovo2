/*
 * read2edge.c
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

#include <stdinc.h>
#include "kmerhash.h"
#include "newhash.h"
#include <extfunc.h>
#include <extvab.h>

static pthread_mutex_t mutex_arc;

long long foundreadcount = 0;
#ifdef MER127
static const Kmer kmerZero = { 0, 0, 0, 0};
#else
static const Kmer kmerZero = { 0, 0 };
#endif

static int buffer_size = 100000000;
//buffer related varibles for chop kmer
static int read_c;
static char ** rcSeq;
static char ** seqBuffer;
static int * lenBuffer;
static char * seqLine;
static int lLineLen = 5000;

// kmer related variables
static int kmer_c;

static Kmer * kmerBuffer;
//static ubyte8 *hashBanBuffer;

static boolean * smallerBuffer;
static int * indexArray;

int nowstep = 0;
int firsttime = 1;

typedef struct fileReadSet
{
	long offset;
	struct fileReadSet * next;
} FILEREADSET;
//record useful reads from parse1readcheck()

static int file_num = 0;
static char ** file_Name;
//record file type:1=fa£¬2=fq
static int * file_type;
//record file max read length
static int * file_maxReadLen;
//record reads offset in file
//static long * offset;

static boolean * nodeBuffer;

FILE * readSeqFile;

static int writeFileNo;
FILE * writeSeqFile;

static void threadRoutine ( void * thrdID );
static void chopKmer4read ( int t, int threadID );
static void thread_wait ( pthread_t * threads );
void Read2edge ( char * libfile, char * graph, int maxk );
void Read2edge2 ( char * libfile, char * graph, int last, int maxk ); //, boolean keepReadFile);
static void parse1readcheck ( int t );
static void searchKmer ( int t, KmerSet2 * kset, int threaID );
static void add1Arc2 ( unsigned int from_ed, unsigned int to_ed, unsigned int weight );
static void sendWorkSignal ( unsigned char SIG, unsigned char * thrdSignals );
static void creatThrds ( pthread_t * threads, PARAMETER * paras );
static void searchKmer1read ( int i, KmerSet2 * kset, int threadID );

struct preArc
{
	unsigned int to_ed;
	unsigned int multiplicity;
	struct preArc * next;
};


struct preArc_array_t
{
	struct preArc ** store_pos;
	unsigned int array_sz;
};
//arc array

struct preArc_array_t arc_arr;
pthread_mutex_t * locks;

/*************************************************
Function:
    put_preArc
Description:
    Puts messages of arc into array.
Input:
    1. arc_arr:     array of arcs
    2. left_id:         from edge id
    3. right_id:        to edge id
    4. added_multi: weight of arc
Output:
    None.
Return:
    None.
*************************************************/
inline void put_preArc ( struct preArc_array_t * arc_arr, unsigned int left_id, unsigned int right_id, unsigned int added_multi )
{
	struct preArc * arc = ( arc_arr->store_pos ) [left_id];

	if ( !arc )
	{
		( arc_arr->store_pos ) [left_id] = ( struct preArc * ) malloc ( sizeof ( struct preArc ) );
		arc = ( arc_arr->store_pos ) [left_id];
		arc->to_ed = right_id;
		arc->multiplicity = added_multi;
		arc->next = NULL;
		return;
	}

	while ( arc )
	{
		if ( arc->to_ed == right_id )
		{
			arc->multiplicity += added_multi;
			return;
		}

		if ( arc->next == NULL ) { break; }

		arc = arc->next;
	}

	arc->next = ( struct preArc * ) malloc ( sizeof ( struct preArc ) );
	arc->next->to_ed = right_id;
	arc->next->multiplicity = added_multi;
	arc->next->next = NULL;
}

//Multi thread.
inline void put_preArc_threaded ( struct preArc_array_t * arc_arr, pthread_mutex_t * locks, unsigned int left_id, unsigned int right_id, unsigned int added_multi )
{
	pthread_mutex_lock ( &locks[left_id] );
	put_preArc ( arc_arr, left_id, right_id, added_multi );
	pthread_mutex_unlock ( &locks[left_id] );
}

//Init.
void init_preArc_array ( struct preArc_array_t * arc_array, unsigned int sz )
{
	arc_array->array_sz = sz;
	arc_array->store_pos = ( struct preArc ** ) calloc ( sz, sizeof ( struct preArc * ) );
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
			{ pthread_join ( threads[i], NULL ); }
}

static void sendWorkSignal ( unsigned char SIG, unsigned char * thrdSignals )
{
	int t;

	for ( t = 0; t < thrd_num; t++ )
		{ thrdSignals[t + 1] = SIG; }

	while ( 1 )
	{
		usleep ( 10 );

		for ( t = 0; t < thrd_num; t++ )
			if ( thrdSignals[t + 1] )
				{ break; }

		if ( t == thrd_num )
			{ break; }
	}
}

static void threadRoutine ( void * para )
{
	PARAMETER * prm;
	int i, t, j, start, finish, k;
	unsigned char id;
	prm = ( PARAMETER * ) para;
	id = prm->threadID;

	//printf("%dth thread with task %d, hash_table %p\n",id,prm.task,prm.hash_table);
	while ( 1 )
	{
		if ( * ( prm->selfSignal ) == 1 )
		{
			for ( i = 0; i < read_c; ++i )
			{
				if ( i % thrd_num != id )
					{ continue; }

				chopKmer4read ( i, id + 1 );
			}

			* ( prm->selfSignal ) = 0;
		}
		else if ( * ( prm->selfSignal ) == 3 )
		{
			* ( prm->selfSignal ) = 0;
			break;
		}
		else if ( * ( prm->selfSignal ) == 5 )
		{
			for ( i = 0; i < read_c; ++i )
			{
				if ( i % thrd_num != id )
					{ continue; }

				searchKmer1read ( i, KmerSetsNew, id );
			}

			* ( prm->selfSignal ) = 0;
		}

		usleep ( 1 );
	}
}

/*************************************************
Function:
    searchKmer1read
Description:
    Searches kmers on reads to add arcs.
Input:
    1. i:               index of reads
    2. kset:            kmer set
    3. threadID:        thread id
Output:
    None.
Return:
    None.
*************************************************/
static void searchKmer1read ( int i, KmerSet2 * kset, int threadID )
{
	kmer_t2 * node1, *node2;
	struct edgeID * edge1, *edge2;
	unsigned int from_ed, to_ed, temp_from, temp_to;
	boolean found;
	ARC * temp = NULL;
	boolean last = false;
	int t;

	for ( t = indexArray[i]; t < indexArray[i + 1] - 1; ++t )
	{
		//get first match
		if ( !last )
			{ found = search_kmerset2 ( kset, kmerBuffer[t], &node1 ); }
		else
		{
			found = true;
			node1 = node2;
		}

		//get next match
		if ( found )
		{
			found = search_kmerset2 ( kset, kmerBuffer[t + step], &node2 );

			if ( found )
			{
				edge1 = node1->edgeId;

				while ( edge1 )
				{
					edge2 = node2->edgeId;

					while ( edge2 )
					{
						temp_from = edge1->edge;
						temp_to = edge2->edge;

						if ( smallerBuffer[t] && smallerBuffer[t + step] )
						{
							if ( ( edge1->flag == 2 || edge1->flag == 4 || edge1->flag == 8 )
							        && ( edge2->flag == 3 || edge2->flag == 4 || edge2->flag == 6 ) )
							{
								{
									if ( temp_from == temp_to || temp_from == getTwinEdge ( temp_to ) )
									{
										pthread_mutex_lock ( &locks[temp_from] );
										edge_array[temp_from].multi = 1;
										pthread_mutex_unlock ( &locks[temp_from] );
										pthread_mutex_lock ( &locks[getTwinEdge ( temp_from )] );
										edge_array[getTwinEdge ( temp_from )].multi = 1;
										pthread_mutex_unlock ( &locks[getTwinEdge ( temp_from )] );
									}
									else
									{
										if ( t == indexArray[i] || t == indexArray[i + 1] - 2 )
										{
											if ( t == indexArray[i] )
												{ put_preArc_threaded ( &arc_arr, locks, temp_from, temp_to, 1 ); }
											else
												{ put_preArc_threaded ( &arc_arr, locks, temp_from, temp_to, 1 ); }
										}
										else
										{
											put_preArc_threaded ( &arc_arr, locks, temp_from, temp_to, 1 );
										}
									}
								}
							}
						}
						else    if ( smallerBuffer[t] && !smallerBuffer[t + step] )
						{
							if ( ( edge1->flag == 2 || edge1->flag == 4 || edge1->flag == 8 )
							        && ( edge2->flag == 1 || edge2->flag == 5 || edge2->flag == 7 ) )
							{
								{
									if ( temp_from == temp_to || temp_from == getTwinEdge ( temp_to ) )
									{
										pthread_mutex_lock ( &locks[temp_from] );
										edge_array[temp_from].multi = 1;
										pthread_mutex_unlock ( &locks[temp_from] );
										pthread_mutex_lock ( &locks[getTwinEdge ( temp_from )] );
										edge_array[getTwinEdge ( temp_from )].multi = 1;
										pthread_mutex_unlock ( &locks[getTwinEdge ( temp_from )] );
									}
									else
									{
										if ( t == indexArray[i] || t == indexArray[i + 1] - 2 )
										{
											if ( t == indexArray[i] )
												{ put_preArc_threaded ( &arc_arr, locks, temp_from, temp_to, 1 ); }
											else
												{ put_preArc_threaded ( &arc_arr, locks, temp_from, temp_to, 1 ); }
										}
										else
										{
											put_preArc_threaded ( &arc_arr, locks, temp_from, temp_to, 1 );
										}
									}
								}
							}
						}
						else    if ( !smallerBuffer[t] && smallerBuffer[t + step] )
						{
							if ( ( edge1->flag == 0 || edge1->flag == 5 || edge1->flag == 9 )
							        && ( edge2->flag == 3 || edge2->flag == 4 || edge2->flag == 6 ) )
							{
								{
									if ( temp_from == temp_to || temp_from == getTwinEdge ( temp_to ) )
									{
										pthread_mutex_lock ( &locks[temp_from] );
										edge_array[temp_from].multi = 1;
										pthread_mutex_unlock ( &locks[temp_from] );
										pthread_mutex_lock ( &locks[getTwinEdge ( temp_from )] );
										edge_array[getTwinEdge ( temp_from )].multi = 1;
										pthread_mutex_unlock ( &locks[getTwinEdge ( temp_from )] );
									}
									else
									{
										if ( t == indexArray[i] || t == indexArray[i + 1] - 2 )
										{
											if ( t == indexArray[i] )
												{ put_preArc_threaded ( &arc_arr, locks, temp_from, temp_to, 1 ); }
											else
												{ put_preArc_threaded ( &arc_arr, locks, temp_from, temp_to, 1 ); }
										}
										else
										{
											put_preArc_threaded ( &arc_arr, locks, temp_from, temp_to, 1 );
										}
									}
								}
							}
						}
						else    if ( !smallerBuffer[t] && !smallerBuffer[t + step] )
						{
							if ( ( edge1->flag == 0 || edge1->flag == 5 || edge1->flag == 9 )
							        && ( edge2->flag == 1 || edge2->flag == 5 || edge2->flag == 7 ) )
							{
								{
									if ( temp_from == temp_to || temp_from == getTwinEdge ( temp_to ) )
									{
										pthread_mutex_lock ( &locks[temp_from] );
										edge_array[temp_from].multi = 1;
										pthread_mutex_unlock ( &locks[temp_from] );
										pthread_mutex_lock ( &locks[getTwinEdge ( temp_from )] );
										edge_array[getTwinEdge ( temp_from )].multi = 1;
										pthread_mutex_unlock ( &locks[getTwinEdge ( temp_from )] );
									}
									else
									{
										if ( t == indexArray[i] || t == indexArray[i + 1] - 2 )
										{
											if ( t == indexArray[i] )
												{ put_preArc_threaded ( &arc_arr, locks, temp_from, temp_to, 1 ); }
											else
												{ put_preArc_threaded ( &arc_arr, locks, temp_from, temp_to, 1 ); }
										}
										else
										{
											put_preArc_threaded ( &arc_arr, locks, temp_from, temp_to, 1 );
										}
									}
								}
							}
						}

						edge2 = edge2->next;
					}

					edge1 = edge1->next;
				}

				if ( firsttime )
				{
					nodeBuffer[t + step] = 1;
				}

				last = true;
			}
			else
			{
				last = false;
			}

			if ( firsttime )
			{
				nodeBuffer[t] = 1;
			}
		}
	}
}

/*************************************************
Function:
    chopKmer4read
Description:
    Chops reads to kmers.
Input:
    1. t:           index of reads
    2. threadID:        thread id
Output:
    None.
Return:
    None.
*************************************************/
static void chopKmer4read ( int t, int threadID )
{
	char * src_seq = seqBuffer[t];
	char * bal_seq = rcSeq[threadID];
	int len_seq = lenBuffer[t];
	int j, bal_j;
	ubyte8 hash_ban, bal_hash_ban;
	Kmer word, bal_word;
	int index;

	if ( len_seq < overlaplen + 1 )
	{
		return;
	}

#ifdef MER127
	word = kmerZero;

	for ( index = 0; index < overlaplen; index++ )
	{
		word = KmerLeftBitMoveBy2 ( word );
		word.low2 |= src_seq[index];
	}

#else
	word = kmerZero;

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
		kmerBuffer[index] = word;
		smallerBuffer[index] = 1;
		index++;
	}
	else
	{
		kmerBuffer[index] = bal_word;
		smallerBuffer[index] = 0;
		index++;
	}

	for ( j = 1; j <= len_seq - overlaplen; j ++ )
	{
		word = nextKmer ( word, src_seq[j - 1 + overlaplen] );
		bal_j = len_seq - j - overlaplen; //  j;
		bal_word = prevKmer ( bal_word, bal_seq[bal_j] );

		if ( KmerSmaller ( word, bal_word ) )
		{
			kmerBuffer[index] = word;
			smallerBuffer[index] = 1;
			index++;
		}
		else
		{
			// complementary node
			kmerBuffer[index] = bal_word;
			smallerBuffer[index] = 0;
			index++;
		}
	}
}

/*************************************************
Function:
    getArcBetween2
Description:
    Get arcs between two edges.
Input:
    1. from_ed:     the from edge
    2. to_ed:       the to edge
Output:
    None.
Return:
    Arc between two edges.
*************************************************/
ARC * getArcBetween2 ( unsigned int from_ed, unsigned int to_ed )
{
	ARC * parc;
	parc = edge_array[from_ed].arcs;

	while ( parc )
	{
		if ( parc->to_ed == to_ed )
			{ return parc; }

		parc = parc->next;
	}

	return parc;
}

/*************************************************
Function:
    add1Arc2
Description:
    Add one arc for two edges.
Input:
    1. from_ed:     the from edge
    2. to_ed:       the to edge
    3. weight:      the weight of arc
Output:
    None.
Return:
    None.
*************************************************/
static void add1Arc2 ( unsigned int from_ed, unsigned int to_ed, unsigned int weight )
{
	unsigned int bal_fe = getTwinEdge ( from_ed );
	unsigned int bal_te = getTwinEdge ( to_ed );

	if ( from_ed > num_ed || to_ed > num_ed || bal_fe > num_ed || bal_te > num_ed )
		{ return; }

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
		{ edge_array[from_ed].arcs->prev = parc; }

	parc->next = edge_array[from_ed].arcs;
	edge_array[from_ed].arcs = parc;

	// A->A'
	if ( bal_te == from_ed )
	{
		parc->bal_arc = parc;
		parc->multiplicity += weight;
		return;
	}

	bal_parc = allocateArc ( bal_fe );
	bal_parc->multiplicity = weight;
	bal_parc->prev = NULL;

	if ( edge_array[bal_te].arcs )
		{ edge_array[bal_te].arcs->prev = bal_parc; }

	bal_parc->next = edge_array[bal_te].arcs;
	edge_array[bal_te].arcs = bal_parc;
	//link them to each other
	parc->bal_arc = bal_parc;
	bal_parc->bal_arc = parc;
}

/*************************************************
Function:
    parse1readcheck
Description:
    Checks whether a read is usablel.
Input:
    1. t:       the index of read
Output:
    None.
Return:
    None.
*************************************************/
static void parse1readcheck ( int t )
{
	unsigned int j;
	unsigned int start, finish;
	boolean found;
	start = indexArray[t];
	finish = indexArray[t + 1];
	boolean readfound = 0;
	nowstep = 1;
	int curr_fileNo;
	FILEREADSET * head, * next, *curr;

	for ( j = start; j < finish - nowstep; ++j )
	{
		found = nodeBuffer[j];

		if ( found )
		{
			found = nodeBuffer[j + nowstep];

			if ( found )
			{
				readfound = 1;
				break;
			}
		}
	}

	if ( readfound )
	{
		++foundreadcount;

		if ( !writeFileNo )
		{
			int index, num_actg;

			for ( index = 0; index < lenBuffer[t]; index++ )
			{
				fprintf ( writeSeqFile, "%c", int2base ( seqBuffer[t][index] ) );
			}

			fprintf ( writeSeqFile, "\n" );
		}
	}
}

//Free.
void free_new()
{
	free_libs();
	free ( file_Name );
	free ( file_type );
	file_type = NULL;
	free ( file_maxReadLen );
	file_maxReadLen = NULL;
}

/*************************************************
Function:
    Read2edge
Description:
    1. Maps reads back to edges.
    2. Outputs selected reads if -r is set.
Input:
    1. libfile:     the reads config file
    2. graph:       the output file prefix
    3. maxk:        the max kmer when using multikmer
Output:
    None.
Return:
    None.
*************************************************/
void Read2edge ( char * libfile, char * graph, int maxk )
{
	long long i;
	char * next_name;
	int maxReadNum, fileNo;
	boolean flag, pairs = 0;
	pthread_t threads[thrd_num];
	unsigned char thrdSignal[thrd_num + 1];
	PARAMETER paras[thrd_num];
	maxReadLen = 0;
	maxNameLen = 256;
	//scan lib info
	scan_libInfo ( libfile );

	if ( !maxReadLen )
		{ maxReadLen = 100; }

	if ( maxk > maxReadLen )
	{
		fprintf ( stderr, "-- Max kmer %d larger than max read length %d, please define a smaller value. --\n", maxk, maxReadLen );
		abort();
	}

	maxReadLen4all = maxReadLen;
	fprintf ( stderr, "In file: %s, max seq len %d, max name len %d.\n",
	          libfile, maxReadLen, maxNameLen );
	int m, n, index;
	file_num = 0;

	for ( m = 0; m < num_libs; m++ )
	{
		if ( lib_array[m].asm_flag == 1 || lib_array[m].asm_flag == 3 )
			file_num += lib_array[m].num_a1_file
			            + lib_array[m].num_a2_file
			            + lib_array[m].num_p_file
			            + lib_array[m].num_q1_file
			            + lib_array[m].num_q2_file
			            + lib_array[m].num_s_a_file
			            + lib_array[m].num_s_q_file;
	}

	file_Name = ( char ** ) ckalloc ( file_num * sizeof ( char * ) );
	file_type = ( int * ) ckalloc ( file_num * sizeof ( int ) );
	file_maxReadLen = ( int * ) ckalloc ( file_num * sizeof ( int ) );
	index = 0;
	//2013-5-14
	int maxReadLenLocal = 0;

	for ( m = 0; m < num_libs; m++ )
	{
		if ( lib_array[m].asm_flag != 1 && lib_array[m].asm_flag != 3 )
			{ continue; }

		//2013-5-14
		if ( lib_array[m].rd_len_cutoff > 0 )
			{ maxReadLenLocal = lib_array[m].rd_len_cutoff < maxReadLen4all ? lib_array[m].rd_len_cutoff : maxReadLen4all; }
		else
			{ maxReadLenLocal = maxReadLen4all; }

		//fa1 fa2
		for ( n = 0; n < lib_array[m].num_a1_file; n++ )
		{
			if ( strlen ( lib_array[m].a1_fname[n] ) > 3 && strcmp ( lib_array[m].a1_fname[n] + strlen ( lib_array[m].a1_fname[n] ) - 3, ".gz" ) == 0 )
			{
				file_type[index] = 3;
				file_maxReadLen[index] = maxReadLenLocal;
				file_Name[index++] = lib_array[m].a1_fname[n];
				file_type[index] = 3;
				file_maxReadLen[index] = maxReadLenLocal;
				file_Name[index++] = lib_array[m].a2_fname[n];
			}
			else
			{
				file_type[index] = 1;
				file_maxReadLen[index] = maxReadLenLocal;
				file_Name[index++] = lib_array[m].a1_fname[n];
				file_type[index] = 1;
				file_maxReadLen[index] = maxReadLenLocal;
				file_Name[index++] = lib_array[m].a2_fname[n];
			}
		}

		//fq1 fq2
		for ( n = 0; n < lib_array[m].num_q1_file; n++ )
		{
			if ( strlen ( lib_array[m].q1_fname[n] ) > 3 && strcmp ( lib_array[m].q1_fname[n] + strlen ( lib_array[m].q1_fname[n] ) - 3, ".gz" ) == 0 )
			{
				file_type[index] = 4;
				file_maxReadLen[index] = maxReadLenLocal;
				file_Name[index++] = lib_array[m].q1_fname[n];
				file_type[index] = 4;
				file_maxReadLen[index] = maxReadLenLocal;
				file_Name[index++] = lib_array[m].q2_fname[n];
			}
			else
			{
				file_type[index] = 2;
				file_maxReadLen[index] = maxReadLenLocal;
				file_Name[index++] = lib_array[m].q1_fname[n];
				file_type[index] = 2;
				file_maxReadLen[index] = maxReadLenLocal;
				file_Name[index++] = lib_array[m].q2_fname[n];
			}
		}

		//fp
		for ( n = 0; n < lib_array[m].num_p_file; n++ )
		{
			if ( strlen ( lib_array[m].p_fname[n] ) > 3 && strcmp ( lib_array[m].p_fname[n] + strlen ( lib_array[m].p_fname[n] ) - 3, ".gz" ) == 0 )
			{
				file_type[index] = 3;
				file_maxReadLen[index] = maxReadLenLocal;
				file_Name[index++] = lib_array[m].p_fname[n];
			}
			else
			{
				file_type[index] = 1;
				file_maxReadLen[index] = maxReadLenLocal;
				file_Name[index++] = lib_array[m].p_fname[n];
			}
		}

		//fa
		for ( n = 0; n < lib_array[m].num_s_a_file; n++ )
		{
			if ( strlen ( lib_array[m].s_a_fname[n] ) > 3 && strcmp ( lib_array[m].s_a_fname[n] + strlen ( lib_array[m].s_a_fname[n] ) - 3, ".gz" ) == 0 )
			{
				file_type[index] = 3;
				file_maxReadLen[index] = maxReadLenLocal;
				file_Name[index++] = lib_array[m].s_a_fname[n];
			}
			else
			{
				file_type[index] = 1;
				file_maxReadLen[index] = maxReadLenLocal;
				file_Name[index++] = lib_array[m].s_a_fname[n];
			}
		}

		//fq
		for ( n = 0; n < lib_array[m].num_s_q_file; n++ )
		{
			if ( strlen ( lib_array[m].s_q_fname[n] ) > 3 && strcmp ( lib_array[m].s_q_fname[n] + strlen ( lib_array[m].s_q_fname[n] ) - 3, ".gz" ) == 0 )
			{
				file_type[index] = 4;
				file_maxReadLen[index] = maxReadLenLocal;
				file_Name[index++] = lib_array[m].s_q_fname[n];
			}
			else
			{
				file_type[index] = 2;
				file_maxReadLen[index] = maxReadLenLocal;
				file_Name[index++] = lib_array[m].s_q_fname[n];
			}
		}
	}

	//init
	next_name = ( char * ) ckalloc ( ( maxNameLen + 1 ) * sizeof ( char ) );
	kmerBuffer = ( Kmer * ) ckalloc ( buffer_size * sizeof ( Kmer ) );
	smallerBuffer = ( boolean * ) ckalloc ( buffer_size * sizeof ( boolean ) );
	pthread_mutex_init ( &mutex_arc, NULL );
	nodeBuffer = ( boolean * ) ckalloc ( buffer_size * sizeof ( boolean ) );
	maxReadNum = buffer_size / ( maxReadLen - overlaplen + 1 );
	seqBuffer = ( char ** ) ckalloc ( maxReadNum * sizeof ( char * ) );
	lenBuffer = ( int * ) ckalloc ( maxReadNum * sizeof ( int ) );
	indexArray = ( int * ) ckalloc ( ( maxReadNum + 1 ) * sizeof ( int ) );
	init_preArc_array ( &arc_arr, num_ed + 2 );
	locks = ( pthread_mutex_t * ) calloc ( arc_arr.array_sz, sizeof ( pthread_mutex_t ) );
	unsigned int ii;

	for ( ii = 0; ii < arc_arr.array_sz; ++ii )
	{
		pthread_mutex_init ( &locks[ii], NULL );
	}

	for ( i = 0; i < maxReadNum; i++ )
		{ seqBuffer[i] = ( char * ) ckalloc ( maxReadLen * sizeof ( char ) ); }

	rcSeq = ( char ** ) ckalloc ( ( thrd_num + 1 ) * sizeof ( char * ) );
	thrdSignal[0] = 0;

	if ( 1 )
	{
		for ( i = 0; i < thrd_num; i++ )
		{
			rcSeq[i + 1] = ( char * ) ckalloc ( maxReadLen * sizeof ( char ) );
			thrdSignal[i + 1] = 0;
			paras[i].threadID = i;
			paras[i].mainSignal = &thrdSignal[0];
			paras[i].selfSignal = &thrdSignal[i + 1];
		}

		creatThrds ( threads, paras );
	}

	if ( 1 )
	{
		rcSeq[0] = ( char * ) ckalloc ( maxReadLen * sizeof ( char ) );
	}

	kmer_c = n_solexa = read_c = i =  readNumBack = gradsCounter = 0;
	fileNo = -1;
	int t0, t1, t2, t3, t4, t5, t6;
	t0 = t1 = t2 = t3 = t4 = t5 = t6 = 0;
	time_t read_start, read_end, time_bef, time_aft;
	time ( &read_start );
	int t;
	char writeSeqName[256];
	writeFileNo = 0;
	sprintf ( writeSeqName, "%s.read", graph );
	writeSeqFile = ckopen ( writeSeqName, "w" );
	int type;
	FILE * file = NULL;
	long pos_seq = 0;

	//parse all reads
	while ( ( flag = read1seqInLibpos ( seqBuffer[read_c], next_name, & ( lenBuffer[read_c] ), //file,
	                                    &fileNo, file_num, file_Name, file_type, file_maxReadLen, &pos_seq ) ) != 0 )
	{
		if ( ( ++i ) % 100000000 == 0 )
			{ fprintf ( stderr, "--- %lldth reads.\n", i ); }

		if ( lenBuffer[read_c] < overlaplen + 1 )
			{ continue; }

		indexArray[read_c] = kmer_c;
		kmer_c += lenBuffer[read_c] - overlaplen + 1;
		read_c++;

		if ( read_c == maxReadNum )
		{
			indexArray[read_c] = kmer_c;
			time ( &read_end );
			t0 += read_end - read_start;
			time ( &time_bef );
			//chop kmer for reads
			sendWorkSignal ( 1, thrdSignal ); //chopKmer4read
			time ( &time_aft );
			t1 += time_aft - time_bef;
			time ( &time_bef );
			//add arc
			sendWorkSignal ( 5, thrdSignal ); //searchKmer1read
			time ( &time_aft );
			t2 += time_aft - time_bef;
			time ( &time_bef );

			//check whether reads is useful
			for ( t = 0; t < read_c; ++t )
			{
				parse1readcheck ( t );
			}

			memset ( nodeBuffer, '\0',  buffer_size * sizeof ( boolean ) );
			time ( &time_aft );
			t3 += time_aft - time_bef;
			kmer_c = 0;
			read_c = 0;
			time ( &read_start );
		}
	}

	//take care of last round
	if ( read_c )
	{
		indexArray[read_c] = kmer_c;
		time ( &read_end );
		t0 += read_end - read_start;
		time ( &time_bef );
		//chop kmer for reads
		sendWorkSignal ( 1, thrdSignal ); //chopKmer4read
		time ( &time_aft );
		t1 += time_aft - time_bef;
		time ( &time_bef );
		//add arc
		sendWorkSignal ( 5, thrdSignal ); //searchKmer1read
		time ( &time_aft );
		t2 += time_aft - time_bef;
		time ( &time_bef );
		struct preArc * parc;
		unsigned int ii;

		for ( ii = 0; ii < arc_arr.array_sz; ++ii )
		{
			parc = ( arc_arr.store_pos ) [ii];

			if ( parc )
			{
				while ( parc )
				{
					add1Arc2 ( ii, parc->to_ed, parc->multiplicity );
					parc = parc->next;
				}
			}
		}

		//check whether reads is useful
		for ( t = 0; t < read_c; ++t )
		{
			parse1readcheck ( t );
		}

		time ( &time_aft );
		t3 += time_aft - time_bef;
	}
	else
	{
		struct preArc * parc;
		unsigned int ii;

		for ( ii = 0; ii < arc_arr.array_sz; ++ii )
		{
			parc = ( arc_arr.store_pos ) [ii];

			if ( parc )
			{
				while ( parc )
				{
					add1Arc2 ( ii, parc->to_ed, parc->multiplicity );
					parc = parc->next;
				}
			}
		}
	}

	fprintf ( stderr, "%lld read(s) processed.\n", i );
	//  fprintf(stderr, "Time spent on reading file: %ds,chop reads: %ds, search kmer: %ds, parse reads: %ds.\n",t0,t1,t2, t3);
	fprintf ( stderr, "Time spent on:\n" );
	fprintf ( stderr, " importing reads: %ds,\n", t0 );
	fprintf ( stderr, " chopping reads to kmers: %ds,\n", t1 );
	fprintf ( stderr, " searching kmers in hash: %ds,\n", t2 );
	fprintf ( stderr, " parsing reads: %ds.\n", t3 );

	if ( foundreadcount )
		{ fprintf ( stderr, "%lld reads available.\n", foundreadcount ); }

	foundreadcount = 0;
	//exit
	sendWorkSignal ( 3, thrdSignal );
	thread_wait ( threads );

	if ( 1 )
	{
		for ( i = 0; i < thrd_num; i++ )
		{
			free ( ( void * ) rcSeq[i + 1] );
		}
	}

	if ( 1 )
	{
		free ( ( void * ) rcSeq[0] );
	}

	free ( ( void * ) rcSeq );
	rcSeq = NULL;

	for ( i = 0; i < maxReadNum; i++ )
		{ free ( ( void * ) seqBuffer[i] ); }

	free ( ( void * ) seqBuffer );
	seqBuffer = NULL;
	free ( ( void * ) lenBuffer );
	lenBuffer = NULL;
	free ( ( void * ) indexArray );
	indexArray = NULL;
	free ( ( void * ) kmerBuffer );
	kmerBuffer = NULL;
	free ( ( void * ) smallerBuffer );
	smallerBuffer = NULL;
	free ( ( void * ) nodeBuffer );
	nodeBuffer = NULL;
	free ( ( void * ) locks );
	struct preArc * temp, *temp_next;

	for ( i = 0; i < arc_arr.array_sz; ++i )
	{
		temp = ( arc_arr.store_pos ) [i];

		while ( temp )
		{
			temp_next = temp->next;
			free ( ( void * ) ( temp ) );
			temp = temp_next;
		}
	}

	free ( ( void * ) arc_arr.store_pos );
	free ( ( void * ) next_name );
	writeFileNo++;
	fclose ( writeSeqFile );
	free_new();
}

int temp_times = 1;

/*************************************************
Function:
    read1seqInNewFile
Description:
    Read one seq in file.
Input:
    1. src_seq:     seq buffer
    2. len_seq:     length of seq
    3. pos_seq:     pos in file
Output:
    None.
Return:
    1 if success.
*************************************************/
static boolean read1seqInNewFile ( char * src_seq, int * len_seq, long * pos_seq )
{
	char c;
	char * str = seqLine;
	int strLen = 0, i;

	/*
	if(temp_times>0)
	{
	    *pos_seq = ftell(readSeqFile);
	}
	*/
	if ( fgets ( str, lLineLen, readSeqFile ) )
	{
		strLen = strlen ( str );
		*len_seq = strLen - 1;

		for ( i = 0; i < strLen - 1; i++ )
		{
			if ( str[i] >= 'a' && str[i] <= 'z' )
			{
				c = base2int ( str[i] - 'a' + 'A' );
				src_seq[i] = c;
			}
			else if ( str[i] >= 'A' && str[i] <= 'Z' )
			{
				c = base2int ( str[i] );
				src_seq[i] = c;
				// after pre-process all the symbles would be a,g,c,t,n in lower or upper case.
			}
			else if ( str[i] == '.' )
			{
				c = base2int ( 'A' );
				src_seq[i] = c;
			}   // after pre-process all the symbles would be a,g,c,t,n in lower or upper case.
		}
	}

	if ( strLen == 0 )
		{ return 0; }
	else
		{ return 1; }
}


/*************************************************
Function:
    Read2edge
Description:
    1. Maps reads back to edges and re-builds arcs between edges.
    2. Outputs selected reads if -r is set.
Input:
    1. libfile:     the reads config file
    2. graph:       the output file prefix
    3. lastTime:        whether it's the last iteration
    4. maxk:        the max kmer when using multikmer
//  5. keepReadFile:        keep tmp reads file that selected for building arcs
Output:
    None.
Return:
    None.
*************************************************/
//void Read2edge2(char *libfile, char *graph, int lastTime, int maxk, boolean keepReadFile)
void Read2edge2 ( char * libfile, char * graph, int lastTime, int maxk )
{
	long long i;
	char * next_name;
	int maxReadNum, fileNo;
	boolean flag, pairs = 0;
	pthread_t threads[thrd_num];
	unsigned char thrdSignal[thrd_num + 1];
	PARAMETER paras[thrd_num];
	maxReadLen = 0;
	maxNameLen = 256;
	//scan lib info
	scan_libInfo ( libfile );

	if ( !maxReadLen )
		{ maxReadLen = 100; }

	if ( maxk > maxReadLen )
	{
		fprintf ( stderr, "-- Max kmer %d larger than max read length %d, please define a smaller value. --\n", maxk, maxReadLen );
		abort();
	}

	maxReadLen4all = maxReadLen;
	fprintf ( stderr, "In file: %s, max seq len %d, max name len %d.\n",
	          libfile, maxReadLen, maxNameLen );
	int m, n, index;
	file_num = 0;

	for ( m = 0; m < num_libs; m++ )
	{
		if ( lib_array[m].asm_flag == 1 || lib_array[m].asm_flag == 3 )
			file_num += lib_array[m].num_a1_file
			            + lib_array[m].num_a2_file
			            + lib_array[m].num_p_file
			            + lib_array[m].num_q1_file
			            + lib_array[m].num_q2_file
			            + lib_array[m].num_s_a_file
			            + lib_array[m].num_s_q_file;
	}

	file_Name = ( char ** ) ckalloc ( file_num * sizeof ( char * ) );
	file_type = ( int * ) ckalloc ( file_num * sizeof ( int ) );
	file_maxReadLen = ( int * ) ckalloc ( file_num * sizeof ( int ) );
	index = 0;
	//2013-5-14
	int maxReadLenLocal = 0;

	for ( m = 0; m < num_libs; m++ )
	{
		if ( lib_array[m].asm_flag != 1 && lib_array[m].asm_flag != 3 )
			{ continue; }

		if ( lib_array[m].rd_len_cutoff > 0 )
			{ maxReadLenLocal = lib_array[m].rd_len_cutoff < maxReadLen4all ? lib_array[m].rd_len_cutoff : maxReadLen4all; }
		else
			{ maxReadLenLocal = maxReadLen4all; }

		//fa1 fa2
		for ( n = 0; n < lib_array[m].num_a1_file; n++ )
		{
			if ( strlen ( lib_array[m].a1_fname[n] ) > 3 && strcmp ( lib_array[m].a1_fname[n] + strlen ( lib_array[m].a1_fname[n] ) - 3, ".gz" ) == 0 )
			{
				file_type[index] = 3;
				file_maxReadLen[index] = maxReadLenLocal;
				file_Name[index++] = lib_array[m].a1_fname[n];
				file_type[index] = 3;
				file_maxReadLen[index] = maxReadLenLocal;
				file_Name[index++] = lib_array[m].a2_fname[n];
			}
			else
			{
				file_type[index] = 1;
				file_maxReadLen[index] = maxReadLenLocal;
				file_Name[index++] = lib_array[m].a1_fname[n];
				file_type[index] = 1;
				file_maxReadLen[index] = maxReadLenLocal;
				file_Name[index++] = lib_array[m].a2_fname[n];
			}
		}

		//fq1 fq2
		for ( n = 0; n < lib_array[m].num_q1_file; n++ )
		{
			if ( strlen ( lib_array[m].q1_fname[n] ) > 3 && strcmp ( lib_array[m].q1_fname[n] + strlen ( lib_array[m].q1_fname[n] ) - 3, ".gz" ) == 0 )
			{
				file_type[index] = 4;
				file_maxReadLen[index] = maxReadLenLocal;
				file_Name[index++] = lib_array[m].q1_fname[n];
				file_type[index] = 4;
				file_maxReadLen[index] = maxReadLenLocal;
				file_Name[index++] = lib_array[m].q2_fname[n];
			}
			else
			{
				file_type[index] = 2;
				file_maxReadLen[index] = maxReadLenLocal;
				file_Name[index++] = lib_array[m].q1_fname[n];
				file_type[index] = 2;
				file_maxReadLen[index] = maxReadLenLocal;
				file_Name[index++] = lib_array[m].q2_fname[n];
			}
		}

		//fp
		for ( n = 0; n < lib_array[m].num_p_file; n++ )
		{
			if ( strlen ( lib_array[m].p_fname[n] ) > 3 && strcmp ( lib_array[m].p_fname[n] + strlen ( lib_array[m].p_fname[n] ) - 3, ".gz" ) == 0 )
			{
				file_type[index] = 3;
				file_maxReadLen[index] = maxReadLenLocal;
				file_Name[index++] = lib_array[m].p_fname[n];
			}
			else
			{
				file_type[index] = 1;
				file_maxReadLen[index] = maxReadLenLocal;
				file_Name[index++] = lib_array[m].p_fname[n];
			}
		}

		//fa
		for ( n = 0; n < lib_array[m].num_s_a_file; n++ )
		{
			if ( strlen ( lib_array[m].s_a_fname[n] ) > 3 && strcmp ( lib_array[m].s_a_fname[n] + strlen ( lib_array[m].s_a_fname[n] ) - 3, ".gz" ) == 0 )
			{
				file_type[index] = 3;
				file_maxReadLen[index] = maxReadLenLocal;
				file_Name[index++] = lib_array[m].s_a_fname[n];
			}
			else
			{
				file_type[index] = 1;
				file_maxReadLen[index] = maxReadLenLocal;
				file_Name[index++] = lib_array[m].s_a_fname[n];
			}
		}

		//fq
		for ( n = 0; n < lib_array[m].num_s_q_file; n++ )
		{
			if ( strlen ( lib_array[m].s_q_fname[n] ) > 3 && strcmp ( lib_array[m].s_q_fname[n] + strlen ( lib_array[m].s_q_fname[n] ) - 3, ".gz" ) == 0 )
			{
				file_type[index] = 4;
				file_maxReadLen[index] = maxReadLenLocal;
				file_Name[index++] = lib_array[m].s_q_fname[n];
			}
			else
			{
				file_type[index] = 2;
				file_maxReadLen[index] = maxReadLenLocal;
				file_Name[index++] = lib_array[m].s_q_fname[n];
			}
		}
	}

	//init
	next_name = ( char * ) ckalloc ( ( maxNameLen + 1 ) * sizeof ( char ) );
	kmerBuffer = ( Kmer * ) ckalloc ( buffer_size * sizeof ( Kmer ) );
	smallerBuffer = ( boolean * ) ckalloc ( buffer_size * sizeof ( boolean ) );
	init_preArc_array ( &arc_arr, num_ed + 2 );
	locks = ( pthread_mutex_t * ) calloc ( arc_arr.array_sz, sizeof ( pthread_mutex_t ) );
	unsigned int ii;

	for ( ii = 0; ii < arc_arr.array_sz; ++ii )
	{
		pthread_mutex_init ( &locks[ii], NULL );
	}

	pthread_mutex_init ( &mutex_arc, NULL );
	firsttime = 0;
	maxReadNum = buffer_size / ( maxReadLen - overlaplen + 1 );
	seqBuffer = ( char ** ) ckalloc ( maxReadNum * sizeof ( char * ) );
	lenBuffer = ( int * ) ckalloc ( maxReadNum * sizeof ( int ) );
	indexArray = ( int * ) ckalloc ( ( maxReadNum + 1 ) * sizeof ( int ) );

	//  offset = (long*)ckalloc((maxReadNum+1)*sizeof(long));

	for ( i = 0; i < maxReadNum; i++ )
		{ seqBuffer[i] = ( char * ) ckalloc ( maxReadLen * sizeof ( char ) ); }

	rcSeq = ( char ** ) ckalloc ( ( thrd_num + 1 ) * sizeof ( char * ) );
	thrdSignal[0] = 0;

	if ( 1 )
	{
		for ( i = 0; i < thrd_num; i++ )
		{
			rcSeq[i + 1] = ( char * ) ckalloc ( maxReadLen * sizeof ( char ) );
			thrdSignal[i + 1] = 0;
			paras[i].threadID = i;
			paras[i].mainSignal = &thrdSignal[0];
			paras[i].selfSignal = &thrdSignal[i + 1];
		}

		creatThrds ( threads, paras );
	}

	if ( 1 )
	{
		rcSeq[0] = ( char * ) ckalloc ( maxReadLen * sizeof ( char ) );
	}

	kmer_c = n_solexa = read_c = i = readNumBack = gradsCounter = 0;
	fileNo = -1;
	int t0, t1, t2, t3, t4, t5, t6;
	t0 = t1 = t2 = t3 = t4 = t5 = t6 = 0;
	time_t read_start, read_end, time_bef, time_aft;
	time ( &read_start );
	//  fprintf(stderr, "Start to read.\n");
	int t;
	char readSeqName[256];
	sprintf ( readSeqName, "%s.read", graph );
	readSeqFile = ckopen ( readSeqName, "r" );
	FILE * file = NULL;
	long pos_seq = 0;
	char tmp[lLineLen];

	if ( maxReadLen > lLineLen )
	{
		lLineLen =  maxReadLen + 1;
		seqLine = ( char * ) ckalloc ( lLineLen * sizeof ( char ) );
	}
	else
	{
		seqLine = tmp;
	}

	//parse all reads
	while ( ( flag = read1seqInNewFile ( seqBuffer[read_c], & ( lenBuffer[read_c] ), &pos_seq ) ) )
	{
		if ( ( ++i ) % 100000000 == 0 )
			{ fprintf ( stderr, "--- %lldth reads.\n", i ); }

		if ( lenBuffer[read_c] < overlaplen + 1 )
			{ continue; }

		//      offset[read_c]=pos_seq;
		indexArray[read_c] = kmer_c;
		kmer_c += lenBuffer[read_c] - overlaplen + 1;
		read_c++;

		if ( read_c == maxReadNum )
		{
			indexArray[read_c] = kmer_c;
			time ( &read_end );
			t0 += read_end - read_start;
			time ( &time_bef );
			//chop kmer
			sendWorkSignal ( 1, thrdSignal ); //chopKmer4read
			time ( &time_aft );
			t1 += time_aft - time_bef;
			time ( &time_bef );
			//add arc
			sendWorkSignal ( 5, thrdSignal ); //searchKmer1read
			time ( &time_aft );
			t2 += time_aft - time_bef;
			time ( &time_bef );
			time ( &time_aft );
			t3 += time_aft - time_bef;
			kmer_c = 0;
			read_c = 0;
			time ( &read_start );
		}
	}

	if ( read_c )
	{
		indexArray[read_c] = kmer_c;
		time ( &read_end );
		t0 += read_end - read_start;
		time ( &time_bef );
		//chop kmer
		sendWorkSignal ( 1, thrdSignal ); //chopKmer4read
		time ( &time_aft );
		t1 += time_aft - time_bef;
		time ( &time_bef );
		//add arc
		sendWorkSignal ( 5, thrdSignal ); //searchKmer1read
		time ( &time_aft );
		t2 += time_aft - time_bef;
		time ( &time_bef );
		struct preArc * parc;
		unsigned int ii;

		for ( ii = 0; ii < arc_arr.array_sz; ++ii )
		{
			parc = ( arc_arr.store_pos ) [ii];

			if ( parc )
			{
				while ( parc )
				{
					add1Arc2 ( ii, parc->to_ed, parc->multiplicity );
					parc = parc->next;
				}
			}
		}

		time ( &time_aft );
		t3 += time_aft - time_bef;
	}
	else
	{
		struct preArc * parc;
		unsigned int ii;

		for ( ii = 0; ii < arc_arr.array_sz; ++ii )
		{
			parc = ( arc_arr.store_pos ) [ii];

			if ( parc )
			{
				while ( parc )
				{
					add1Arc2 ( ii, parc->to_ed, parc->multiplicity );
					parc = parc->next;
				}
			}
		}
	}

	fprintf ( stderr, "%lld read(s) processed.\n", i );
	//  fprintf(stderr, "Time spent on reading file: %ds,chop reads: %ds, search kmer: %ds, parse reads: %ds.\n",t0,t1,t2, t3);
	fprintf ( stderr, "Time spent on:\n" );
	fprintf ( stderr, " importing reads: %ds,\n", t0 );
	fprintf ( stderr, " chopping reads to kmers: %ds,\n", t1 );
	fprintf ( stderr, " searching kmers in hash: %ds,\n", t2 );
	fprintf ( stderr, " parsing reads: %ds.\n", t3 );

	if ( foundreadcount )
		{ fprintf ( stderr, "%lld reads available.\n", foundreadcount ); }

	foundreadcount = 0;
	//exit
	sendWorkSignal ( 3, thrdSignal );
	thread_wait ( threads );

	if ( 1 )
	{
		for ( i = 0; i < thrd_num; i++ )
		{
			free ( ( void * ) rcSeq[i + 1] );
		}
	}

	if ( 1 )
	{
		free ( ( void * ) rcSeq[0] );
	}

	free ( ( void * ) rcSeq );
	rcSeq = NULL;

	for ( i = 0; i < maxReadNum; i++ )
		{ free ( ( void * ) seqBuffer[i] ); }

	free ( ( void * ) seqBuffer );
	seqBuffer = NULL;
	free ( ( void * ) lenBuffer );
	lenBuffer = NULL;
	free ( ( void * ) indexArray );
	indexArray = NULL;
	free ( ( void * ) kmerBuffer );
	kmerBuffer = NULL;
	free ( ( void * ) smallerBuffer );
	smallerBuffer = NULL;
	//  free((void*)offset);
	//  offset=NULL;
	free ( ( void * ) locks );
	struct preArc * temp, *temp_next;

	for ( i = 0; i < arc_arr.array_sz; ++i )
	{
		temp = ( arc_arr.store_pos ) [i];

		while ( temp )
		{
			temp_next = temp->next;
			free ( ( void * ) ( temp ) );
			temp = temp_next;
		}
	}

	free ( ( void * ) arc_arr.store_pos );
	free ( ( void * ) next_name );
	fclose ( readSeqFile );

	if ( !lastTime )
	{
		//      if(!keepReadFile)
		remove ( readSeqName );
	}

	free_new();

	if ( maxReadLen > lLineLen )
	{
		free ( ( void * ) seqLine );
	}
}


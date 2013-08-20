/*
 * build_preArc.cpp
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
#include "build_preArc.h"
#include "seq_util.h"
#include "global.h"
#include "multi_threads.h"

#include <math.h>
#include "bam.h"
#include "faidx.h"
#include "knetfile.h"
#include "sam_view.h"
#include "xcurses.h"
#include "zlib.h"
#include "bgzf.h"
#include "glf.h"
#include "kstring.h"
#include "razf.h"
#include "sam_header.h"
#include "zconf.h"

#include "sam.h"

#include "sparse_kmer.h"

//programming guide:

//step1:read edge & build  vertex hash
//step1_include:
//a.create id for every edge & it's rev_comp id
//b.build vertex for every edge (done...)

//step2:read reads & chop kmers

//step3:search vertex with the kmers & mark the found position & get the Vertex
//step3_extend:if found two or more vetex in one read , record the path to solve tiny repeat

//step4:compare the left kmer right kmer to vertex's edge_kmer

//step5: if left & right both match build a preArc or increase the multiplicity by 1



void init_vertex_hash ( vertex_hash2 * v_ht, size_t sz )
{
	v_ht->ht_sz = sz;
	v_ht->store_pos = ( vertex2 ** ) calloc ( sz, sizeof ( vertex2 * ) );
}


/*************************************************
Function:
    build_vertexes
Description:
    1. Reads sequence from *.sparse.edge.
    2. Builds vertexes by cutting the edge sequence's end kmers.
Input:
    1. v_ht:        vertex hashtable
    2. K_size:      kmer size
    3. edge_file:       edge file
Output:
    None.
Return:
    None.
*************************************************/
void build_vertexes ( vertex_hash2 * v_ht, int K_size, char * edge_file )
{
	FILE * fp;
	kmer_t2 from_kmer, to_kmer;
	size_t line_len, edge_len_left;
	int edge_len;
	int cvg;
	bool bal_ed;//回文为0
	const int BUFF_LEN = 1024;
	char line[BUFF_LEN];
	char str[32];
	char to_buff[BUFF_LEN];//buffer 2k edge seq  BUFF_LEN>4*K_size
	int processed = 0; //0表示未处理 1表示 处理完毕 2 表示已经处理了from_vertex ,to_vertex 还没有处理
	size_t edge_id = 0;
	fp = fopen ( edge_file, "r" );

	if ( !fp )
	{
		fprintf ( stderr, "ERROR: Cannot open edge_file %s. Now exit to system...\n", edge_file );
		exit ( -1 );
	}

	vertex2 * v_tmp;
	edge_starter2 * e_tmp;
	int is_found;
	bool is_left;

	while ( fgets ( line, BUFF_LEN, fp ) != NULL )
	{
		//debug<<"processed "<<processed<<endl;
		if ( line[0] == '>' ) //get one edge length, from vertex, to vertex,cvg,bal
		{
			if ( processed == 1 && bal_ed )
			{
				edge_id++;//如果不是回文
			}

#ifdef _63MER_
			sscanf ( line + 7, "%d,%llx %llx ,%llx %llx ,%s %d,%d", &edge_len,
			         & ( from_kmer.kmer ) [0], & ( from_kmer.kmer ) [1], & ( to_kmer.kmer ) [0], & ( to_kmer.kmer ) [1], str, &cvg, &bal_ed ); // from_kmer to_kmer is of no use here
#endif
#ifdef _127MER_
			sscanf ( line + 7, "%d,%llx %llx %llx %llx ,%llx %llx %llx %llx ,%s %d,%d", &edge_len,
			         & ( from_kmer.kmer ) [0], & ( from_kmer.kmer ) [1], & ( from_kmer.kmer ) [2], & ( from_kmer.kmer ) [3],
			         & ( to_kmer.kmer ) [0], & ( to_kmer.kmer ) [1], & ( to_kmer.kmer ) [2], & ( to_kmer.kmer ) [3], str, &cvg, &bal_ed ); // from_kmer to_kmer is of no use here
#endif
			edge_len_left = K_size + edge_len;
			processed = 0;
			edge_id++;// current edge positive strand id
			//debug<<line<<"edge_id "<<edge_id<<endl;
		}
		else
		{
			if ( processed == 0 )
			{
				line_len = strlen ( line );

				if ( line[line_len - 1] == '\n' )
				{
					line[line_len - 1] = '\0';
					line_len --;
				}

				if ( edge_len_left - line_len == 0 ) //edge completely loaded
				{
					//do all process
					process_edge ( v_ht, K_size, line, line_len, 1, edge_id, bal_ed );
					processed = 1;
					edge_len_left = 0;
					continue;
				}
				else  //edge partly loaded at the first time.
				{
					if ( line_len < 2 * K_size ) //line_len < 2*K_size &&edge_len_left - line_len > 0
					{
						fprintf ( stderr, "ERROR:it won't happen in 63mer/127mer\n" );
						exit ( 1 );
					}
					else
					{
						process_edge ( v_ht, K_size, line, line_len, 2, edge_id, bal_ed );
						processed = 2;
						edge_len_left -= line_len;

						if ( edge_len_left >= 2 * K_size )
						{
							//no need to buf the to kmer seq
						}
						else if ( edge_len_left < 2 * K_size )
						{
							//to_buff[100];/ copy the last 2K char of line  to to_buff  already no '\n'
							strcpy ( to_buff, line + ( line_len - 2 * K_size ) );
						}
						else
						{
							fprintf ( stderr, "ERROR: in cal the edge_len_left!!\n" );
							exit ( 1 );
						}
					}
				}
			}
			else if ( processed == 2 )
			{
				//if(line[0]=='\n') continue;
				line_len = strlen ( line );

				if ( line[line_len - 1] == '\n' )
				{
					line[line_len - 1] = '\0';
					line_len --;
				}

				edge_len_left -= line_len;

				if ( edge_len_left == 0 ) //load the complete edge sequence
				{
					//process the to kmer
					if ( line_len >= 2 * K_size )
					{
						process_edge ( v_ht, K_size, line, line_len, 3, edge_id, bal_ed );
					}
					else
					{
						//need to use the to_buff
						int buf_len = strlen ( to_buff );
						strcpy ( to_buff + buf_len, line );
						buf_len = strlen ( to_buff );
						process_edge ( v_ht, K_size, to_buff, buf_len, 3, edge_id, bal_ed );
					}

					processed = 1;
					continue;
				}
				else
				{
					if ( edge_len_left >= 2 * K_size )
					{
						//no need to buf the to kmer seq
					}
					else if ( edge_len_left < 2 * K_size )
					{
						//to_buff[100];/ copy the last 2K char of line  to to_buff
						strcpy ( to_buff, line + ( line_len - 2 * K_size ) );
					}
					else
					{
						fprintf ( stderr, "ERROR: in cal the edge_len_left!!\n" );
						exit ( 1 );
					}
				}
			}
			else
			{
				if ( line[0] == '\n' ) { continue; } //当len = 1023时

				fprintf ( stderr, "ERROR: in cal the status_processed !! %d \n", processed );
				exit ( 1 );
			}
		}
	}

	fclose ( fp );
}


/*************************************************
Function:
    process_edge
Description:
    It builds vetexes  from one or part of one edge sequence.
Input:
    1. v_ht:        hashtable
    2. K_size:      kmer size
    3. seq:     edge sequence
    4. len:     edge length
    5. type:        1: process head and tail;  2: process head ; 3:process tail
    6. edge_id:     edge id
    7. bal_edge:        0:palindrome 1:else
Output:
    None.
Return:
    None.
*************************************************/
static void process_edge ( vertex_hash2 * v_ht, int K_size, char * seq, int len, int type, size_t edge_id, bool bal_edge )
{
	kmer_t2 vertex_kmer;
	kmer_t2 edge_kmer;
	vertex2 * v_tmp;
	edge_starter2 * e_tmp;
	int is_found;
	bool is_left;
	int edge_kmer_len;

	switch ( type )
	{
		case 1: //process all ..
			//process the head
			get_kmer_from_seq ( seq, len, K_size, 0, &vertex_kmer );

			if ( len <= K_size + gap ) //get the last kmer
			{
				get_kmer_from_seq ( seq, len, K_size, len - K_size, &edge_kmer );
				edge_kmer_len = len - K_size;
			}
			else
			{
				//get_kmer_from_seq(seq, len, K_size, K_size,&edge_kmer);
				get_kmer_from_seq ( seq, len, K_size, gap, &edge_kmer );
				edge_kmer_len = gap;
			}

			is_left = 0;//right
			v_tmp = put_vertex ( v_ht, vertex_kmer, is_found );
			put_edge ( v_tmp, edge_kmer,  is_left, edge_kmer_len, edge_id );
			reverseCompKmer ( &vertex_kmer, K_size );
			reverseCompKmer ( &edge_kmer, K_size );
			is_left = 1;//left
			v_tmp = put_vertex ( v_ht, vertex_kmer, is_found );
			put_edge ( v_tmp, edge_kmer,  is_left, edge_kmer_len, edge_id + bal_edge );
			//process the tail
			get_kmer_from_seq ( seq, len, K_size, len - K_size, &vertex_kmer );

			if ( len <= K_size + gap ) //get the first kmer
			{
				get_kmer_from_seq ( seq, len, K_size, 0, &edge_kmer );
				edge_kmer_len = len - K_size;
			}
			else
			{
				get_kmer_from_seq ( seq, len, K_size, len - K_size - gap, &edge_kmer );
				edge_kmer_len = gap;
			}

			is_left = 1;
			v_tmp = put_vertex ( v_ht, vertex_kmer, is_found );
			put_edge ( v_tmp, edge_kmer,  is_left, edge_kmer_len, edge_id );
			reverseCompKmer ( &vertex_kmer, K_size );
			reverseCompKmer ( &edge_kmer, K_size );
			is_left = 0;//right
			v_tmp = put_vertex ( v_ht, vertex_kmer, is_found );
			put_edge ( v_tmp, edge_kmer,  is_left, edge_kmer_len, edge_id + bal_edge );
			break;
		case 2:
			//process only the  head
			get_kmer_from_seq ( seq, len, K_size, 0, &vertex_kmer );

			if ( len <= K_size + gap )
			{
				get_kmer_from_seq ( seq, len, K_size, len - K_size, &edge_kmer );
				edge_kmer_len = len - K_size;
			}
			else
			{
				get_kmer_from_seq ( seq, len, K_size, gap, &edge_kmer );
				edge_kmer_len = gap;
			}

			is_left = 0;//right
			v_tmp = put_vertex ( v_ht, vertex_kmer, is_found );
			put_edge ( v_tmp, edge_kmer,  is_left, edge_kmer_len, edge_id );
			reverseCompKmer ( &vertex_kmer, K_size );
			reverseCompKmer ( &edge_kmer, K_size );
			is_left = 1;//left
			v_tmp = put_vertex ( v_ht, vertex_kmer, is_found );
			put_edge ( v_tmp, edge_kmer,  is_left, edge_kmer_len, edge_id + bal_edge );
			break;
		case 3:
			//process only the tail
			get_kmer_from_seq ( seq, len, K_size, len - K_size, &vertex_kmer );

			if ( len <= K_size + gap )
			{
				get_kmer_from_seq ( seq, len, K_size, 0, &edge_kmer );
				edge_kmer_len = len - K_size;
			}
			else
			{
				get_kmer_from_seq ( seq, len, K_size, len - K_size - gap, &edge_kmer );
				edge_kmer_len = gap;
			}

			is_left = 1;
			v_tmp = put_vertex ( v_ht, vertex_kmer, is_found );
			put_edge ( v_tmp, edge_kmer,  is_left, edge_kmer_len, edge_id );
			reverseCompKmer ( &vertex_kmer, K_size );
			reverseCompKmer ( &edge_kmer, K_size );
			is_left = 0;//right
			v_tmp = put_vertex ( v_ht, vertex_kmer, is_found );
			put_edge ( v_tmp, edge_kmer,  is_left, edge_kmer_len, edge_id + bal_edge );
			break;
		default:
			fprintf ( stderr, "ERROR: wrong process type in process_edge()\n" );
			exit ( 1 );
	}
}


static vertex2 * put_vertex ( vertex_hash2 * v_ht, kmer_t2 vertex_kmer, int & is_found ) //63 127 differ fixed
{
	uint64_t hv = MurmurHash64A ( vertex_kmer.kmer, sizeof ( kmer_t2 ), 0 ); //hash value
	uint64_t idx = ( size_t ) ( hv % v_ht->ht_sz );
	vertex2 * ver = ( v_ht->store_pos ) [idx];

	if ( !ver )
	{
		( v_ht->store_pos ) [idx] = ( vertex2 * ) malloc ( sizeof ( vertex2 ) );
		ver = ( v_ht->store_pos ) [idx];
		ver->kmer_t2 = vertex_kmer;
		ver->left = NULL;
		ver->right = NULL;
		ver->next = NULL;
		is_found = 0;
		return ver;
	}

	while ( ver )
	{
		if ( kmerCompare ( & ( ver->kmer_t2 ), &vertex_kmer ) == 0 )
		{
			is_found = 1;
			return ver;
		}

		if ( ver->next == NULL ) { break; }

		ver = ver->next;
	}

	is_found = 0;
	ver->next = ( vertex2 * ) malloc ( sizeof ( vertex2 ) );
	ver->next->kmer_t2 = vertex_kmer;
	ver->next->left = NULL;
	ver->next->right = NULL;
	ver->next->next = NULL;
	return ver->next;
}

static void put_edge ( vertex2 * ver, kmer_t2 edge_kmer, bool is_left, int len, size_t edge_id ) //fixed
{
	edge_starter2 * tmp = NULL;

	if ( is_left )
	{
		if ( !ver->left )
		{
			ver->left = ( edge_starter2 * ) malloc ( sizeof ( edge_starter2 ) );
			ver->left->edge_kmer = edge_kmer;
			ver->left->edge_id = edge_id;
			ver->left->len = len;//record the length of edge (1~k)
			ver->left->next = NULL;
			return;
		}

		tmp = ver->left;
	}
	else
	{
		if ( !ver->right )
		{
			ver->right = ( edge_starter2 * ) malloc ( sizeof ( edge_starter2 ) );
			ver->right->edge_kmer = edge_kmer;
			ver->right->edge_id = edge_id;
			ver->right->len = len;//record the length of edge (1~k)
			ver->right->next = NULL;
			return;
		}

		tmp = ver->right;
	}

	while ( tmp->next ) //because there are no two edges equal attached with one node ...
	{
		tmp = tmp->next;
	}

	tmp->next = ( edge_starter2 * ) malloc ( sizeof ( edge_starter2 ) );
	tmp->next->edge_kmer = edge_kmer;
	tmp->next->edge_id = edge_id;
	tmp->next->len = len;//record the length of edge (1~k)
	tmp->next->next = NULL;
}




static vertex2 * search_vertex ( vertex_hash2 * v_ht, kmer_t2 * vertex_kmer ) //fixed ...
{
	uint64_t hv = MurmurHash64A ( vertex_kmer->kmer, sizeof ( kmer_t2 ), 0 ); //hash value
	uint64_t idx = ( size_t ) ( hv % v_ht->ht_sz );
	vertex2 * ver = ( v_ht->store_pos ) [idx];

	while ( ver )
	{
		if ( kmerCompare ( & ( ver->kmer_t2 ), vertex_kmer ) == 0 )
		{
			return ver;
		}

		ver = ver->next;
	}

	return NULL;
}



void init_preArc_array ( preArc_array * arc_array, size_t sz ) //63 127 same
{
	arc_array->array_sz = sz;
	arc_array->store_pos = ( preArc ** ) calloc ( sz, sizeof ( preArc * ) );
}


static void chop_kmers ( const char * read, int len, int K_size, kmer_t2 * kmer_array, int kmer_array_len, int & kmer_num )
{
	if ( len <= K_size )
	{
		kmer_num = 0;
		return ;
	}

	kmer_num = len - K_size + 1;

	if ( kmer_num > kmer_array_len )
	{
		fprintf ( stderr, "ERROR: the kmer_array_len is not enough! %d\n", kmer_num );
		exit ( 1 );
	}

	kmer_t2 kmer;

	for ( int i = 0; i < kmer_num; ++i )             //optimize later
	{
		get_kmer_from_seq ( read, len, K_size, i, &kmer );
		kmer_array[i] = kmer;
	}
}

/*************************************************
Function:
    put_preArc_threaded
Description:
    Stores one preArc into preArc_array
Input:
    1. arc_arr:     preArc array
    2. locks:       the locks array for modifying preArc array
    3. left_id:     left edge id
    4. right_id:        right edge id
    5. added_multi:     added weigth
Output:
    None.
Return:
    None.
*************************************************/
static inline void put_preArc_threaded ( preArc_array * arc_arr, pthread_spinlock_t * locks, size_t left_id, size_t right_id, int added_multi )
{
	pthread_spin_lock ( &locks[left_id] );
	put_preArc ( arc_arr, left_id, right_id, added_multi );
	pthread_spin_unlock ( &locks[left_id] );
}


static inline void put_preArc ( preArc_array * arc_arr, size_t left_id, size_t right_id, int added_multi )
{
	preArc * arc = ( arc_arr->store_pos ) [left_id];

	if ( !arc )
	{
		( arc_arr->store_pos ) [left_id] = ( preArc * ) malloc ( sizeof ( preArc ) );
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

	arc->next = ( preArc * ) malloc ( sizeof ( preArc ) );
	arc->next->to_ed = right_id;
	arc->next->multiplicity = added_multi;
	arc->next->next = NULL;
}


static inline void put_preArc ( preArc_array * arc_arr, size_t left_id, size_t right_id )
{
	preArc * arc = ( arc_arr->store_pos ) [left_id];

	if ( !arc )
	{
		( arc_arr->store_pos ) [left_id] = ( preArc * ) malloc ( sizeof ( preArc ) );
		arc = ( arc_arr->store_pos ) [left_id];
		arc->to_ed = right_id;
		arc->multiplicity = 1;
		arc->next = NULL;
		return;
	}

	while ( arc )
	{
		if ( arc->to_ed == right_id )
		{
			arc->multiplicity++;
			return;
		}

		if ( arc->next == NULL ) { break; }

		arc = arc->next;
	}

	arc->next = ( preArc * ) malloc ( sizeof ( preArc ) );
	arc->next->to_ed = right_id;
	arc->next->multiplicity = 1;
	arc->next->next = NULL;
}


void output_preArcs ( preArc_array * arc_arr, char * outfile )
{
	FILE * fp;
	fp = fopen ( outfile, "w" );

	if ( !fp )
	{
		fprintf ( stderr, "ERROR: can't create file %s in output_preArc\n", outfile );
		exit ( 1 );
	}

	preArc * parc;

	for ( size_t i = 0; i < arc_arr->array_sz; ++i )
	{
		parc = ( arc_arr->store_pos ) [i];

		if ( parc )
		{
			fprintf ( fp, "%u", i );
			int j = 0;

			while ( parc )
			{
				j++;
				fprintf ( fp, " %u %u", parc->to_ed, parc->multiplicity );
				parc = parc->next;

				if ( parc && j % 4 == 0 )
				{
					fprintf ( fp, "\n" );
					fprintf ( fp, "%u", i );
				}
			}

			fprintf ( fp, "\n" );
		}
	}

	fclose ( fp );
}


static void free_vertex ( vertex2 * tmp )
{
	edge_starter2 * edge_s, *edge_s2;
	edge_s = tmp->left;

	while ( edge_s )
	{
		edge_s2 = edge_s;
		edge_s = edge_s->next;
		free ( edge_s2 );
	}

	edge_s = tmp->right;

	while ( edge_s )
	{
		edge_s2 = edge_s;
		edge_s = edge_s->next;
		free ( edge_s2 );
	}

	free ( tmp );
}

void free_vertex_hash ( vertex_hash2 * v_ht )
{
	vertex2 * tmp, *tmp2;

	for ( size_t i = 0; i < v_ht->ht_sz; ++i )
	{
		tmp = ( v_ht->store_pos ) [i];

		while ( tmp )
		{
			tmp2 = tmp;
			tmp = tmp->next;
			free_vertex ( tmp2 );
		}
	}

	free ( v_ht->store_pos );
}

/*************************************************
Function:
    process_1read_preArc
Description:
    This is the core function for building preArcs.
    1. Chops one read into kmers.
    2. Searches the kmers in vertex hash.
    3. Aligns the vertex's kmer-edge sequences to the read sequence on both sides.
    4. Constructs preArcs according the mapping result on both sides of a vertex.

    @since r53:
    5. add -R support, solves tiny repeat.
Input:
    1. arc_arr:     preArc array
    2. locks:       locks array
    3. v_ht:        vertex hash
    4. K_size:      kmer size
    5. cut_off_len:     cut off length
    6. read:        read
Output:
    None.
Return:
    None.
*************************************************/
void process_1read_preArc ( preArc_array * arc_arr, pthread_spinlock_t * locks, int thread_id, vertex_hash2 * v_ht, int K_size, int cut_off_len, const char * read )
{
	const int BUFF_LEN = 1024;
	kmer_t2 kmers[BUFF_LEN];
	int kmer_array_len = cut_off_len - K_size + 1;
	int kmer_num ;
	vertex2 * v_tmp;
	edge_starter2 * e_tmp;
	size_t left_id;
	size_t right_id;
	int left_found = 0, right_found = 0;
	int edge_len;
	//update
	//int map_len;
	//int shortest_maplen = 0;
	//add for -R solving tiny repeats
	unsigned int path[128];
	unsigned int counter = 0;
	//int read_len,i=0;
	int read_len = strlen ( read );
	/*
	while(read[i]!='\0'){
	    i++;
	}
	read_len = i;
	//read_len = strlen(read);
	if(read[read_len-1]=='\n'){
	    read[read_len-1]='\0';
	    read_len--;
	}*/

	if ( read_len > cut_off_len ) { read_len = cut_off_len; }

	kmer_array_len = read_len - K_size + 1;
	chop_kmers ( read, read_len, K_size, kmers, kmer_array_len, kmer_num );

	for ( int i = 1; i < kmer_num - 1; ++i ) //search every kmer exclude the begin and end kmer
	{
		v_tmp = search_vertex ( v_ht, &kmers[i] );

		if ( v_tmp ) //found
		{
			//search left edge kmer got left id
			e_tmp = v_tmp->left;

			while ( e_tmp )
			{
				edge_len = e_tmp->len;

				if ( edge_len <= i )
				{
					if ( kmerCompare ( & ( kmers[i - edge_len] ), & ( e_tmp->edge_kmer ) ) == 0 )
					{
						left_id = e_tmp->edge_id;

						if ( left_found )
						{
							fprintf ( stderr, "ERROR: left edge id found already !new found id %llu \n", left_id );
							fprintf ( stderr, "i:%d ,edge_len:%d\n", i, edge_len );
							printKmerSeq ( & ( kmers[i - edge_len] ), K_size, stderr );
							printKmerSeq ( & ( e_tmp->edge_kmer ), K_size, stderr );
							exit ( 1 );
						};

						left_found = 1;

						break;
					}
				}
				else
				{
					kmer_t2 read_edge = kmers[0];

					if ( K_size > i )
					{
						kmerMoveRight ( &read_edge, K_size - i );
					}

					kmer_t2 KMER_FILTER;
					initKmerFilter ( i, &KMER_FILTER );
					kmer_t2 edge_kmer = e_tmp->edge_kmer;

					if ( K_size > edge_len )
					{
						kmerMoveRight ( &edge_kmer, K_size - edge_len );
					}

					kmerAnd ( &read_edge, &KMER_FILTER );
					kmerAnd ( &edge_kmer, &KMER_FILTER );

					if ( kmerCompare ( &read_edge, &edge_kmer ) == 0 )
					{
						left_found++;
						left_id = e_tmp->edge_id;

						if ( left_found == 2 )
						{
							//debug_build<<"can't distinct which left edge\n";
							break;
						}
					}
				}

				e_tmp = e_tmp->next;
			}

			//update maplen_control
			/*
			if(edge_len >= shortest_maplen){
			    if(map_len < shortest_maplen) left_found = 0;
			}else{
			    if(map_len != edge_len) left_found = 0;
			}*/

			if ( left_found != 1 ) {left_found = 0; right_found = 0; continue;} //not found or multi found

			//todo : aln  if  left_found = 0  ... find the best
			//search right edge kmer got right id
			e_tmp = v_tmp->right;

			while ( e_tmp )
			{
				edge_len = e_tmp->len;

				if ( edge_len <= kmer_num - 1 - i )
				{
					if ( kmerCompare ( & ( kmers[i + edge_len] ), & ( e_tmp->edge_kmer ) ) == 0 )
					{
						right_id = e_tmp->edge_id;

						if ( right_found )
						{
							fprintf ( stderr, "ERROR: right edge id found already, new found id %llu !\n", right_id );
							fprintf ( stderr, "i:%d ,edge_len:%d\n", i, edge_len );
							printKmerSeq ( & ( kmers[i + edge_len] ), K_size, stderr );
							printKmerSeq ( & ( e_tmp->edge_kmer ), K_size, stderr );
							exit ( 1 );
						};

						right_found = 1;

						break;
					}
				}
				else
				{
					int read_edge_len = ( kmer_num - 1 - i );
					kmer_t2 KMER_FILTER;
					initKmerFilter ( read_edge_len, &KMER_FILTER );
					kmer_t2 read_edge = kmers[kmer_num - 1];
					kmerAnd ( &read_edge, &KMER_FILTER );
					kmer_t2 edge_kmer = e_tmp->edge_kmer;

					if ( edge_len > read_edge_len )
					{
						kmerMoveRight ( &edge_kmer, ( edge_len - read_edge_len ) );
					}

					kmerAnd ( &edge_kmer, &KMER_FILTER );

					if ( kmerCompare ( &read_edge, &edge_kmer ) == 0 )
					{
						right_found++;
						right_id = e_tmp->edge_id;

						if ( right_found == 2 )
						{
							//debug_build<<"can't distinct which right edge\n";
							break;
						}
					}
				}

				e_tmp = e_tmp->next;
			}

			//update map_len control
			/*
			if(edge_len >= shortest_maplen){
			    if(map_len < shortest_maplen) right_found = 0;
			}else{
			    if(map_len != edge_len) right_found = 0;
			}*/

			if ( right_found != 1 ) {left_found = 0; right_found = 0; continue;}

			//todo : aln  if  right_found = 0  ... find the best
			//if(left_found == 1 && right_found ==1)
			//store this preArc
			//preArc_array *arc_arr
			put_preArc_threaded ( arc_arr, locks, left_id, right_id, 1 );

			//constructing the path ...
			if ( solve )
			{
				if ( counter == 0 )
				{
					counter = 2;
					path[1] = left_id;
					path[2] = right_id;
				}
				else if ( counter <= 100 )
				{
					if ( path[counter] == left_id )
					{
						path[++counter] = right_id;
					}
					else
					{
						path[++counter] = left_id;
						path[++counter] = right_id;
					}
				}
			}

			//end ...
			left_found = 0;
			right_found = 0;
		}
	}

	//add to path buffer , if full filled ,output it
	if ( solve )
	{
		if ( counter >= 3 && counter <= 100 )
		{
			path[0] = counter;
			int tmp = is_full ( path_buffer[thread_id] );

			if ( tmp == 1 )
			{
				//output it
				output_edge_path_buffer_locked ( path_buffer[thread_id], path_fp, &file_lock );
			}
			else if ( tmp == -1 )
			{
				//error status
				fprintf ( stderr, "ERROR: path buffer overflow!! system exit .\n" );
				exit ( -1 );
			}

			put_path_2_buffer ( path_buffer[thread_id],  path );
		}
	}
}


void free_preArc_array ( preArc_array * arc_array )
{
	preArc * tmp, *tmp2;

	for ( size_t i = 0; i < arc_array->array_sz; ++i )
	{
		tmp = ( arc_array->store_pos ) [i];

		while ( tmp )
		{
			tmp2 = tmp;
			tmp = tmp->next;
			free ( tmp2 );
		}
	}

	free ( arc_array->store_pos );
}



/*************************************************
Function:
    build_preArc_threaded
Description:
    This is the main entry for building preArcs.
Input:
    1. arc_arr:     preArc array
    2. v_ht:        vertex hash
    3. K_size:      kmer size
    4. cut_off_len:     cut off length
    5. in_filenames_vt:     input reads file names
    6. thread_num:      thread number
Output:
    None.
Return:
    None.
*************************************************/
void build_preArc_threaded ( preArc_array * arc_arr, vertex_hash2 * v_ht, int K_size, int cut_off_len, vector<string> *in_filenames_vt, int thread_num )
{
	//create main io thread
	int read_buf_sz = 102400 * thrd_num_s;
	read_buf0 = new string[read_buf_sz];
	read_buf1 = new string[read_buf_sz];
	io_stat0 = 1; //must be one, if io_stat0 =0 ,the io thread will work immediately
	io_stat1 = 1;
	io_ready = 0;
	io_para_main io_para_mains;
	io_para_mains.read_buf_sz = read_buf_sz;
	io_para_mains.in_filenames_vt = in_filenames_vt;
	pthread_t io_thread;
	int temp;

	//fprintf(stderr,"Creating main io thread ...\n");
	if ( ( temp = pthread_create ( &io_thread, NULL, run_io_thread_main, &io_para_mains ) ) != 0 )
	{
		fprintf ( stderr, "ERROR: failed creating main io thread.\n" );
		exit ( -1 );
	}

	fprintf ( stderr, "1 io thread initialized.\n" );
	//create work threads ..
	//fprintf(stderr,"Creating work threads ...\n");
	pthread_t threads[thrd_num_s];
	unsigned char thrdSignal[thrd_num_s + 1];
	PARAMETER paras[thrd_num_s];
	locks = ( pthread_spinlock_t * ) calloc ( arc_arr->array_sz, sizeof ( pthread_spinlock_t ) );

	//init as unlock stat ..
	for ( size_t i = 0; i < arc_arr->array_sz; ++i )
	{
		locks[i] = 1;
	}

	for ( int k = 0; k < thrd_num_s; k++ )
	{
		thrdSignal[k + 1] = 0;
		paras[k].threadID = k;
		paras[k].mainSignal = &thrdSignal[0];
		paras[k].selfSignal = &thrdSignal[k + 1];
		paras[k].ht = NULL;
		paras[k].preArcs = arc_arr;
		paras[k].v_ht = v_ht;
		paras[k].cut_off_len = cut_off_len;
		paras[k].K_size = K_size;
		paras[k].gap = gap;
	}

	creatThrds ( threads, paras );
	thrdSignal[0] = 0;

	//run it
	while ( 1 )
	{
		sendIOWorkSignal();

		while ( io_ready == 0 ) {usleep ( 1 );}

		if ( io_ready )
		{
			sendWorkSignal ( 12, thrdSignal );
		}

		if ( io_ready == 2 )
		{
			//fprintf(stderr,"All reads have been processed!\n");
			break;
		}
	}

	sendWorkSignal ( 3, thrdSignal );
	thread_wait ( threads );
	delete [] read_buf0;
	delete [] read_buf1;
	free ( ( void * ) locks );
	free_vertex_hash ( v_ht );
}


/*************************************************
Function:
    create_edge_path_buffer
Description:
    Creates an edge_path_buffer struct dynamicly.
Input:
    None.
Output:
    None.
Return:
    an edge_path_buffer pointer to  heap
*************************************************/
edge_path_buffer * create_edge_path_buffer
( unsigned int * mark_on_edge,
  pthread_spinlock_t * locks,
  unsigned long long buff_size,
  unsigned int max_path_length )
{
	if ( ! ( mark_on_edge && locks ) )
	{
		fprintf ( stderr, "ERROR: The initial mark_on_edge array or locks are not valid! Exit System ...\n" );
		exit ( -1 );
	}

	edge_path_buffer * new_buffer = ( edge_path_buffer * ) calloc ( 1, sizeof ( edge_path_buffer ) );
	new_buffer->mark_on_edge = mark_on_edge;
	new_buffer->locks = locks;
	new_buffer->buff_size = buff_size;
	new_buffer->max_path_length = max_path_length;
	new_buffer->filled_num = 0;
	new_buffer->path_buffer = NULL;
	unsigned int ** tmp;
	tmp = ( unsigned int ** ) calloc ( buff_size, sizeof ( unsigned int * ) );

	for ( size_t i = 0; i < buff_size; i++ )
	{
		tmp[i] = ( unsigned int * ) calloc ( max_path_length, sizeof ( unsigned int ) );
	}

	new_buffer->path_buffer = tmp;
	return new_buffer;
}

/*************************************************
Function:
    create_edge_path_buffer
Description:
    free the path_buffer in an edge_path_buffer struct
Input:
    None
Output:
    None
Return:
    None
*************************************************/
void destory_edge_path_buffer ( struct edge_path_buffer * buffer )
{
	unsigned int ** tmp = buffer->path_buffer;

	for ( size_t i = 0; i < buffer->buff_size; i++ )
	{
		free ( ( void * ) ( tmp[i] ) );
	}

	free ( ( void * ) tmp );
	buffer->filled_num = 0;
}

void clear_edge_path_buffer ( struct edge_path_buffer * buffer )
{
	unsigned int ** tmp = buffer->path_buffer;

	for ( size_t i = 0; i < buffer->buff_size; i++ )
	{
		memset ( tmp[i], 0, buffer->max_path_length * sizeof ( unsigned int ) );
	}

	buffer->filled_num = 0;
}


void output_edge_path_buffer ( struct edge_path_buffer * buffer, FILE * path_file )
{
	if ( debug )
	{
		static size_t times = 0, total = 0;
		total += buffer->filled_num;
		fprintf ( stderr, "call output_edge_path_buffer %lu %lu\n", times++, total );
	}

	if ( !path_file )
	{
		fprintf ( stderr, "ERROR: The path_file is not avilable!\n" );
		exit ( -1 );
	}

	unsigned int counter;
	unsigned int ** tmp = buffer->path_buffer;

	for ( size_t i = 0; i < buffer->filled_num; i++ )
	{
		counter = tmp[i][0];
		fwrite ( &counter, sizeof ( char ), 1, path_file );
		fwrite ( tmp[i] + 1, sizeof ( unsigned int ), ( int ) counter, path_file );
	}

	buffer->filled_num = 0;
}


/*************************************************
Function:
    output_edge_path_buffer_locked
Description:
    Output the buffer to a file for multi-threading.
Input:
    None.
Output:
    None.
Return:
    None.
*************************************************/
void output_edge_path_buffer_locked ( struct edge_path_buffer * buffer, FILE * path_file, pthread_mutex_t * file_mutex )
{
	static size_t times = 0, total = 0;;

	if ( !path_file )
	{
		fprintf ( stderr, "ERROR: The path_file is not avilable!\n" );
		exit ( -1 );
	}

	unsigned int counter;
	unsigned int ** tmp = buffer->path_buffer;
	pthread_mutex_lock ( file_mutex );

	if ( debug )
	{
		total += buffer->filled_num;
		fprintf ( stderr, "call output_edge_path_buffer_locked %lu %lu\n", times++, total );
	}

	for ( size_t i = 0; i < buffer->filled_num; i++ )
	{
		counter = tmp[i][0];
		fwrite ( &counter, sizeof ( char ), 1, path_file );
		fwrite ( tmp[i] + 1, sizeof ( unsigned int ), ( int ) counter, path_file );
	}

	pthread_mutex_unlock ( file_mutex );
	buffer->filled_num = 0;
}



/*************************************************
Function:
    put_path_2_buffer
Description:
    Stores the path to buffer and update mark_on_edge array at the same time.
Input:
    None.
Output:
    None.
Return:
    None.
*************************************************/
int put_path_2_buffer ( struct edge_path_buffer * buffer, unsigned int * path )
{
	if ( debug )
	{
		static size_t times = 0;
		static pthread_spinlock_t lock = 1;
		pthread_spin_lock ( &lock );
		fprintf ( stderr, "call put_path_2_buffer %lu\n", times++ );
		pthread_spin_unlock ( &lock );
	}

	unsigned long long pos = buffer->filled_num;

	if ( pos >= buffer->buff_size )
	{
		return -1;
	}

	memcpy ( ( buffer->path_buffer ) [pos], path, buffer->max_path_length * sizeof ( unsigned int ) );

	for ( unsigned int i = 1; i < path[0]; i++ )
	{
		pthread_spin_lock ( ( buffer->locks ) + path[i] );
		( ( buffer->mark_on_edge ) [path[i]] ) ++;
		pthread_spin_unlock ( ( buffer->locks ) + path[i] );
	}

	buffer->filled_num++;
	return 1;
}



int is_full ( struct edge_path_buffer * buffer )
{
	if ( buffer->filled_num == buffer->buff_size )
	{
		return 1;
	}
	else if ( buffer->filled_num < buffer->buff_size )
	{
		return 0;
	}
	else
	{
		return -1;
	}
}


void clear_status ( struct edge_path_buffer * buffer )
{
	buffer->filled_num = 0;
}
















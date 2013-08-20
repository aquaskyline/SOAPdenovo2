/*
 * inc/build_preArc.h
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
#ifndef _BUILD_PREARC_H
#define _BUILD_PREARC_H

#include "core.h"
#include "global.h"
#include "stdinc.h"

struct edge_starter2
{
	struct kmer_t2 edge_kmer;
	uint64_t edge_id: 32, len: 6; //make sure that left always be end & right always be start
	edge_starter2 * next;

};
struct vertex2
{
	struct kmer_t2 kmer_t2;
	edge_starter2 * left;
	edge_starter2 * right;
	vertex2 * next;
};
struct vertex_hash2
{
	struct vertex2 ** store_pos;
	size_t ht_sz;
};

struct preArc
{
	unsigned int to_ed;
	unsigned int multiplicity;
	struct preArc * next;
};

struct preArc_array
{
	struct preArc ** store_pos;
	size_t array_sz;
};

//public methods
void init_vertex_hash ( vertex_hash2 * v_ht, size_t sz );
void build_vertexes ( vertex_hash2 * v_ht, int K_size, char * edge_file );
void free_vertex_hash ( vertex_hash2 * v_ht );


void init_preArc_array ( preArc_array * arc_array, size_t sz );
void build_preArc_threaded ( preArc_array * arc_arr, vertex_hash2 * v_ht, int K_size, int cut_off_len, vector<string> *in_filenames_vt, int thread_num );
void output_preArcs ( preArc_array * arc_arr, char * outfile );
void free_preArc_array ( preArc_array * arc_array );



//local structs ...
struct io_para
{
	//char **buf0;
	//char **buf1;
	int * io_stat0;
	int * io_stat1;
	int * read_num0;
	int * read_num1;

	int * finished_arr0;
	int * finished_arr1;

	//FILE *fp;

	vector<string> *in_filenames_vt;

	int read_buf_sz;
	int read_buf_len;

	int thread_num;
};


struct process_para
{
	//char **buf0;
	//char **buf1;
	int * io_stat0;
	int * io_stat1;
	int * read_num0;
	int * read_num1;

	int * finished_arr0;
	int * finished_arr1;

	preArc_array * preArcs;//change preArc**  to  preArc* for spin_lock version
	pthread_spinlock_t * locks; //...

	int thread_id;
	int thread_num;

	vertex_hash2 * v_ht;
	int K_size;
	int cut_off_len;

};

void process_1read_preArc ( preArc_array * arc_arr, pthread_spinlock_t * locks, int thread_id, vertex_hash2 * v_ht, int K_size, int cut_off_len, const char * read );


//static methods
static void process_edge ( vertex_hash2 * v_ht, int K_size, char * seq, int len, int type, size_t edge_id, bool bal_edge );
static vertex2 * put_vertex ( vertex_hash2 * v_ht, kmer_t2 vertex_kmer, int & is_found );
static void put_edge ( vertex2 * ver, kmer_t2 edge_kmer, bool is_left, int len, size_t edge_id );
static vertex2 * search_vertex ( vertex_hash2 * v_ht, kmer_t2 * vertex_kmer );
static void free_vertex ( vertex2 * tmp );

static void  get_kmer ( const char * seq, int len, int K_size, int pos, kmer_t2 & kmer );
static void chop_kmers ( char * read, int len, int K_size, kmer_t2 * kmer_array, int kmer_array_len, int & kmer_num );

static void * run_io_thread ( void * arg );
static void * run_process_thread ( void * arg );



static void put_preArc ( preArc_array * arc_arr, size_t left_id, size_t right_id, int added_multi );
static void put_preArc_threaded ( preArc_array * arc_arr, pthread_spinlock_t * locks, size_t left_id, size_t right_id, int added_multi );

//add for solving repeat

struct edge_path_buffer
{
	unsigned int * mark_on_edge; //The mark on edge array, record the times of occurrence for each edge  and it's revers complement
	pthread_spinlock_t * locks; //the locks for multi threads access and modification to mark_on_edge
	unsigned int ** path_buffer; //buffered the paths for out put, (the first unsigned int is the length of the path.)
	unsigned int max_path_length; // the max length for each path   //set to 255 default ...
	unsigned long long buff_size;      // the  max path number the buffer can sotre
	unsigned long long filled_num;    // the filled number of the buffer
};


struct edge_path_buffer * create_edge_path_buffer ( unsigned int * mark_on_edge, pthread_spinlock_t * locks, unsigned long long buff_size, unsigned int max_path_length );
void destory_edge_path_buffer ( struct edge_path_buffer * buffer );
void clear_edge_path_buffer ( struct edge_path_buffer * buffer );
void output_edge_path_buffer ( struct edge_path_buffer * buffer, FILE * path_file );
void output_edge_path_buffer_locked ( struct edge_path_buffer * buffer, FILE * path_file, pthread_mutex_t * file_mutex );


int put_path_2_buffer ( struct edge_path_buffer * buffer, unsigned int * path );


void clear_status ( struct edge_path_buffer * buffer );
int is_full ( struct edge_path_buffer * buffer );



#endif







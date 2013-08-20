/*
 * global.cpp
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

#include "global.h"

char shortrdsfile[256];//in lib file
char graphfile[256];//out prefix

int NodeCovTh = 1;
int EdgeCovTh = 1; //BE :build edge ,BPA:build preArc

int K_size = 23;
int gap = 15;
uint64_t GenomeSize = 0;
int solve = 0; //default not solving repeats
int run_mode = 0;

int thrd_num_s = 8;

size_t * edge_cnt_total = NULL; //used int lock strategy
size_t * bucket_count_total = NULL; //used in lock strategy


//for io thread @see io_func.h
//for io thread
string * seq_t = NULL;
int io_ready;       //0 ready to work 1 working 2 *seq_t ready 3 end reading signal
int read_num = 0;   //the read num in *seq_t

string * read_buf0 = NULL;
string * read_buf1 = NULL;
int io_stat0 = 1; //must be one, if io_stat0 =0 ,the io thread will work immediately
int io_stat1 = 1;

size_t reads_all_num = 0;

int max_rd_len = 0;
int min_rd_len = 100000;


//for the hashing lock strategy ...
pthread_spinlock_t * locks;



// solving tiny repeats, temporarily using global vars to implements this feature
unsigned int * mark_on_edge = NULL;
pthread_spinlock_t * s_locks = NULL;
struct edge_path_buffer ** path_buffer = NULL;
unsigned long long buff_size = 1024;
unsigned int max_path_length = 128; //max_path_length-1 is the real max path length, because the first int of buffer record the path length

FILE * mark_fp = NULL; //
FILE * path_fp = NULL; //
pthread_mutex_t file_lock;//

int debug = 0 ;



















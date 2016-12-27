/*
 * inc/global.h
 *
 * Copyright (c) 2008-2016 Ruibang Luo <aquaskyline.com>.
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

#ifndef _GLOBAL_H
#define _GLOBAL_H

#include "stdinc.h"
#include "io_func.h"

#ifdef __APPLE__
#include "spinLock.h"
#endif


//the definitions @see global.cpp

extern char shortrdsfile[256];//in lib file
extern char graphfile[256];//out prefix

extern int NodeCovTh;
extern int EdgeCovTh;//BE :build edge ,BPA:build preArc

extern int K_size;
extern int gap;
extern uint64_t GenomeSize;
extern int solve; // solve reapeat or not
extern int run_mode;

extern int thrd_num_s;

extern size_t *edge_cnt_total ;  //used int lock strategy
extern size_t *bucket_count_total ;  //used in lock strategy


//for io thread @see io_func.h
extern string *seq_t;
extern int io_ready;
extern int read_num;

extern string *read_buf0;
extern string *read_buf1;
extern int io_stat0; //must be one, if io_stat0 =0 ,the io thread will work immediately
extern int io_stat1;

extern size_t reads_all_num;

extern int max_rd_len;
extern int min_rd_len;

//for the lock strategy ...
extern pthread_spinlock_t *locks;


extern unsigned int *mark_on_edge;
extern pthread_spinlock_t *s_locks;
extern struct edge_path_buffer **path_buffer;
extern unsigned long long buff_size;
extern unsigned int max_path_length; //max_path_length-1 is the real max path length, because the first int of buffer record the path length

extern FILE *mark_fp ;  //
extern FILE *path_fp ;  //
extern pthread_mutex_t file_lock;//

extern int debug;

#endif





































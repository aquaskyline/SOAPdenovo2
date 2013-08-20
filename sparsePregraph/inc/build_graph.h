/*
 * inc/build_graph.h
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


#ifndef __BUILD_GRAPH_H
#define __BUILD_GRAPH_H


#include "stdinc.h"
#include "core.h"
//for called

void run_process_threaded ( struct hashtable2 * ht, pthread_spinlock_t * locks, int K_size, int gap, size_t read_num, int thrd_num, int thrd_id, int round );
void SwitchBuckets ( hashtable2 * ht2, int K_size );
void SavingSparseKmerGraph2 ( hashtable2 * ht, char * outfile );
void LoadingSparseKmerGraph2 ( hashtable2 * ht, char * outfile );

#endif





































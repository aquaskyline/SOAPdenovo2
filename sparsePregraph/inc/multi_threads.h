/*
 * inc/multi_threads.h
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

#ifndef _MULTI_THREADS_H
#define _MULTI_THREADS_H
#include "stdinc.h"

typedef struct parameter
{
	unsigned char threadID;
	struct hashtable2 * ht;
	struct preArc_array * preArcs; //for building preArc ...
	struct vertex_hash2 * v_ht; //for building preArc ...
	int cut_off_len;
	int K_size;
	int gap;
	unsigned char * mainSignal;
	unsigned char * selfSignal;
} PARAMETER;

void creatThrds ( pthread_t * threads, PARAMETER * paras );
void * threadRoutine ( void * para );
void thread_wait ( pthread_t * threads );
void sendWorkSignal ( unsigned char SIG, unsigned char * thrdSignals );

#endif
























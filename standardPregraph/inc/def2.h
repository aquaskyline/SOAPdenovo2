/*
 * inc/def2.h
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

#ifndef _DEF2
#define _DEF2
typedef char boolean;
typedef long long IDnum;
typedef double Time;
typedef long long Coordinate;
// Fibonacci heaps used mainly in Tour Bus
typedef struct fibheap FibHeap;
typedef struct fibheap_el FibHeapNode;
typedef struct dfibheap DFibHeap;
typedef struct dfibheap_el DFibHeapNode;
//Memory manager
typedef struct block_start
{
	struct block_start * next;
} BLOCK_START;

typedef struct recycle_mark
{
	struct recycle_mark * next;
} RECYCLE_MARK;

typedef struct mem_manager
{
	BLOCK_START * block_list;
	int index_in_block;
	int items_per_block;
	size_t item_size;
	RECYCLE_MARK * recycle_list;
	unsigned long long counter;
} MEM_MANAGER;

struct dfibheap_el
{
	int dfhe_degree;
	boolean dfhe_mark;
	DFibHeapNode * dfhe_p;
	DFibHeapNode * dfhe_child;
	DFibHeapNode * dfhe_left;
	DFibHeapNode * dfhe_right;
	Time dfhe_key;
	unsigned int dfhe_data;//void *dfhe_data;
};
#endif

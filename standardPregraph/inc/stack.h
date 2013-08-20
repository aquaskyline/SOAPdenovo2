/*
 * inc/stack.h
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

#ifndef __STACK__
#define __STACK__

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

typedef struct block_starter
{
	struct block_starter * prev;
	struct block_starter * next;
} BLOCK_STARTER;

typedef struct stack
{
	BLOCK_STARTER * block_list;
	int index_in_block;
	int items_per_block;
	int item_c;
	size_t item_size;
	BLOCK_STARTER * block_backup;
	int index_backup;
	int item_c_backup;
} STACK;

void stackBackup ( STACK * astack );
void stackRecover ( STACK * astack );
void * stackPush ( STACK * astack );
void * stackPop ( STACK * astack );
void freeStack ( STACK * astack );
void emptyStack ( STACK * astack );
STACK * createStack ( int num_items, size_t unit_size );


#endif

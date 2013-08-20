/*
 * mem_manager.c
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

MEM_MANAGER * createMem_manager ( int num_items, size_t unit_size )
{
	MEM_MANAGER * mem_Manager = ( MEM_MANAGER * ) ckalloc ( 1 * sizeof ( MEM_MANAGER ) );
	mem_Manager->block_list = NULL;
	mem_Manager->items_per_block = num_items;
	mem_Manager->item_size = unit_size;
	mem_Manager->recycle_list = NULL;
	mem_Manager->counter = 0;
	return mem_Manager;
}

void freeMem_manager ( MEM_MANAGER * mem_Manager )
{
	BLOCK_START * ite_block, *temp_block;

	if ( !mem_Manager )
	{
		return;
	}

	ite_block = mem_Manager->block_list;

	while ( ite_block )
	{
		temp_block = ite_block;
		ite_block = ite_block->next;
		free ( ( void * ) temp_block );
	}

	free ( ( void * ) mem_Manager );
}

void * getItem ( MEM_MANAGER * mem_Manager )
{
	RECYCLE_MARK * mark; //this is the type of return value
	BLOCK_START * block;

	if ( !mem_Manager )
	{
		return NULL;
	}

	if ( mem_Manager->recycle_list )
	{
		mark = mem_Manager->recycle_list;
		mem_Manager->recycle_list = mark->next;
		return mark;
	}

	mem_Manager->counter++;

	if ( !mem_Manager->block_list || mem_Manager->index_in_block == mem_Manager->items_per_block )
	{
		//pthread_mutex_lock(&gmutex);
		block = ckalloc ( sizeof ( BLOCK_START ) + mem_Manager->items_per_block * mem_Manager->item_size );
		//mem_Manager->counter += sizeof(BLOCK_START)+mem_Manager->items_per_block*mem_Manager->item_size;
		//pthread_mutex_unlock(&gmutex);
		block->next = mem_Manager->block_list;
		mem_Manager->block_list = block;
		mem_Manager->index_in_block = 1;
		return ( RECYCLE_MARK * ) ( ( void * ) block + sizeof ( BLOCK_START ) );
	}

	block = mem_Manager->block_list;
	return ( RECYCLE_MARK * ) ( ( void * ) block + sizeof ( BLOCK_START ) + mem_Manager->item_size * ( mem_Manager->index_in_block++ ) );
}

void returnItem ( MEM_MANAGER * mem_Manager, void * item )
{
	RECYCLE_MARK * mark;
	mark = item;
	mark->next = mem_Manager->recycle_list;
	mem_Manager->recycle_list = mark;
}

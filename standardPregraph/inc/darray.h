/*
 * inc/darray.h
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

#ifndef __DARRAY__
#define __DARRAY__

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

typedef struct dynamic_array
{
	void * array;
	long long array_size;
	size_t item_size;
	long long item_c;
} DARRAY;

void * darrayPut ( DARRAY * darray, long long index );
void * darrayGet ( DARRAY * darray, long long index );
DARRAY * createDarray ( int num_items, size_t unit_size );
void freeDarray ( DARRAY * darray );
void emptyDarray ( DARRAY * darray );

#endif


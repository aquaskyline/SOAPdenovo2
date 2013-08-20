/*
 * dfibHeap.c
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

#include <stdlib.h>
#include <stdio.h>

#include  "def2.h"
#include "dfib.h"

// Return number of elements stored in heap
IDnum getDFibHeapSize ( DFibHeap * heap )
{
	return dfibheap_getSize ( heap );
}

// Constructor
// Memory allocated
DFibHeap * newDFibHeap ()
{
	return dfh_makekeyheap ();
}

// Add new node into heap with a key, and a pointer to the specified node
DFibHeapNode * insertNodeIntoDHeap ( DFibHeap * heap, Time key, unsigned int node )
{
	DFibHeapNode * res;
	res = dfh_insertkey ( heap, key, node );
	return res;
}

// Replaces the key for a given node
Time replaceKeyInDHeap ( DFibHeap * heap, DFibHeapNode * node, Time newKey )
{
	Time res;
	res = dfh_replacekey ( heap, node, newKey );
	return res;
}


/*************************************************
Function:
    removeNextNodeFromDHeap
Description:
    Removes the edge from DHeap, then returns it.
Input:
    1. heap :       the heap
Output:
    None.
Return:
    The key.
*************************************************/
unsigned int removeNextNodeFromDHeap ( DFibHeap * heap )
{
	unsigned int node;
	node = ( unsigned int ) dfh_extractmin ( heap );
	return node;
}

// Destructor
void destroyDHeap ( DFibHeap * heap )
{
	dfh_deleteheap ( heap );
}

// Replace the node pointed to by a heap node
void replaceValueInDHeap ( DFibHeapNode * node, unsigned int newValue )
{
	dfh_replacedata ( node, newValue );
}

// Remove unwanted node
void destroyNodeInDHeap ( DFibHeapNode * node, DFibHeap * heap )
{
	dfh_delete ( heap, node );
}

Time getKey ( DFibHeapNode * node )
{
	return dfibheap_el_getKey ( node );
}

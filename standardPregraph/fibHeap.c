/*
 * fibHeap.c
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

#include "fib.h"

// Constructor
// Memory allocated
FibHeap * newFibHeap ()
{
	return fh_makekeyheap ();
}

// Add new node into heap with a key, and a pointer to the specified node
FibHeapNode * insertNodeIntoHeap ( FibHeap * heap, Coordinate key, unsigned int node )
{
	return fh_insertkey ( heap, key, node );
}

// Returns smallest key in heap
Coordinate minKeyOfHeap ( FibHeap * heap )
{
	return fh_minkey ( heap );
}

// Replaces the key for a given node
Coordinate replaceKeyInHeap ( FibHeap * heap, FibHeapNode * node, Coordinate newKey )
{
	return fh_replacekey ( heap, node, newKey );
}

// Removes the node with the shortest key, then returns it.
unsigned int removeNextNodeFromHeap ( FibHeap * heap )
{
	return ( unsigned int ) fh_extractmin ( heap );
}

boolean IsHeapEmpty ( FibHeap * heap )
{
	return fh_isempty ( heap );
}

// Destructor
void destroyHeap ( FibHeap * heap )
{
	fh_deleteheap ( heap );
}

// Replace the node pointed to by a heap node
void replaceValueInHeap ( FibHeapNode * node, unsigned int newValue )
{
	fh_replacedata ( node, newValue );
}

// Remove unwanted node
void destroyNodeInHeap ( FibHeapNode * node, FibHeap * heap )
{
	fh_delete ( heap, node );
}

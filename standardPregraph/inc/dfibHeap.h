/*
 * inc/dfibHeap.h
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

#ifndef _DFIBHEAP_H_
#define _DFIBHEAP_H_

DFibHeap * newDFibHeap();

DFibHeapNode * insertNodeIntoDHeap ( DFibHeap * heap, Time key, unsigned int node );

Time replaceKeyInDHeap ( DFibHeap * heap, DFibHeapNode * node, Time newKey );

unsigned int removeNextNodeFromDHeap ( DFibHeap * heap );

void destroyDHeap ( DFibHeap * heap );

boolean HasMin ( DFibHeap * h );

void replaceValueInDHeap ( DFibHeapNode * node, unsigned int newValue );

void * destroyNodeInDHeap ( DFibHeapNode * node, DFibHeap * heap );

IDnum getDFibHeapSize ( DFibHeap * heap );

Time getKey ( DFibHeapNode * node );
#endif

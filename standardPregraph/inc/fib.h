/*
 * inc/fib.h
 *
 * This file is part of SOAPdenovo.
 *
 */

/*-
 * Copyright 1997, 1998-2003 John-Mark Gurney.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 */

#ifndef _FIB_H_
#define _FIB_H_

//#include "globals.h"
#include <stdio.h>
#include "def2.h"

typedef Coordinate ( *voidcmp ) ( unsigned int , unsigned int );

/* functions for key heaps */
boolean fh_isempty ( FibHeap * );
FibHeap * fh_makekeyheap ( void );
FibHeapNode * fh_insertkey ( FibHeap *, Coordinate, unsigned int );
Coordinate fh_minkey ( FibHeap * );
Coordinate fh_replacekey ( FibHeap *, FibHeapNode *, Coordinate );
unsigned int fh_replacekeydata ( FibHeap *, FibHeapNode *, Coordinate, unsigned int );

/* functions for unsigned int * heaps */
FibHeap * fh_makeheap ( void );
voidcmp fh_setcmp ( FibHeap *, voidcmp );
unsigned int fh_setneginf ( FibHeap *, unsigned int );
FibHeapNode * fh_insert ( FibHeap *, unsigned int );

/* shared functions */
unsigned int fh_extractmin ( FibHeap * );
unsigned int fh_min ( FibHeap * );
unsigned int fh_replacedata ( FibHeapNode *, unsigned int );
unsigned int fh_delete ( FibHeap *, FibHeapNode * );
void fh_deleteheap ( FibHeap * );
FibHeap * fh_union ( FibHeap *, FibHeap * );

#endif              /* _FIB_H_ */

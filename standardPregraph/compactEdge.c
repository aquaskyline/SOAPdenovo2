/*
 * compactEdge.c
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

void copyEdge ( unsigned int source, unsigned int target )
{
	edge_array[target].from_vt = edge_array[source].from_vt;
	edge_array[target].to_vt = edge_array[source].to_vt;
	edge_array[target].length = edge_array[source].length;
	edge_array[target].cvg = edge_array[source].cvg;
	edge_array[target].multi = edge_array[source].multi;
	edge_array[target].flag = edge_array[source].flag;

	if ( edge_array[target].seq )
	{
		free ( ( void * ) edge_array[target].seq );
	}

	edge_array[target].seq = edge_array[source].seq;
	edge_array[source].seq = NULL;
	edge_array[target].arcs = edge_array[source].arcs;
	edge_array[source].arcs = NULL;
	edge_array[target].markers = edge_array[source].markers;
	edge_array[source].markers = NULL;
	edge_array[target].deleted = edge_array[source].deleted;
}

//move edge from source to target
void edgeMove ( unsigned int source, unsigned int target )
{
	unsigned int bal_source, bal_target;
	ARC * arc;
	copyEdge ( source, target );
	bal_source = getTwinEdge ( source );

	//bal_edge
	if ( bal_source != source )
	{
		bal_target = target + 1;
		copyEdge ( bal_source, bal_target );
		edge_array[target].bal_edge = 2;
		edge_array[bal_target].bal_edge = 0;
	}
	else
	{
		edge_array[target].bal_edge = 1;
		bal_target = target;
	}

	//take care of the arcs
	arc = edge_array[target].arcs;

	while ( arc )
	{
		arc->bal_arc->to_ed = bal_target;
		arc = arc->next;
	}

	if ( bal_target == target )
	{
		return;
	}

	arc = edge_array[bal_target].arcs;

	while ( arc )
	{
		arc->bal_arc->to_ed = target;
		arc = arc->next;
	}
}

/*************************************************
Function:
    compactEdgeArray
Description:
    Compacts the edge array by removing deleted edges.
Input:
    None.
Output:
    None.
Return:
    None.
*************************************************/
void compactEdgeArray ()
{
	unsigned int i;
	unsigned int validCounter = 0;
	unsigned int bal_ed;
	fprintf ( stderr, "Before compacting, %d edge(s) existed.\n", num_ed );

	for ( i = 1; i <= num_ed; i++ )
	{
		if ( edge_array[i].deleted )
		{
			continue;
		}

		validCounter++;

		if ( i == validCounter )
		{
			continue;
		}

		bal_ed = getTwinEdge ( i );
		edgeMove ( i, validCounter );

		if ( bal_ed != i )
		{
			i++;
			validCounter++;
		}
	}

	num_ed = validCounter;
	fprintf ( stderr, "After compacting, %d edge(s) left.\n", num_ed );
}

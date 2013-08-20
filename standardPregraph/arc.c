/*
 * arc.c
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

#define preARCBLOCKSIZE 100000

void createPreArcMemManager ()
{
	prearc_mem_manager = createMem_manager ( preARCBLOCKSIZE, sizeof ( preARC ) );
}

void prlDestroyPreArcMem ()
{
	if ( !preArc_mem_managers )
	{
		return;
	}

	int i;

	for ( i = 0; i < thrd_num; i++ )
	{
		freeMem_manager ( preArc_mem_managers[i] );
	}

	free ( ( void * ) preArc_mem_managers );
	preArc_mem_managers = NULL;
}

void destroyPreArcMem ()
{
	freeMem_manager ( prearc_mem_manager );
	prearc_mem_manager = NULL;
}

preARC * prlAllocatePreArc ( unsigned int edgeid, MEM_MANAGER * manager )
{
	preARC * newArc;
	newArc = ( preARC * ) getItem ( manager );
	newArc->to_ed = edgeid;
	newArc->multiplicity = 1;
	newArc->next = NULL;
	return newArc;
}

preARC * allocatePreArc ( unsigned int edgeid )
{
	arcCounter++;
	preARC * newArc;
	newArc = ( preARC * ) getItem ( prearc_mem_manager );
	newArc->to_ed = edgeid;
	newArc->multiplicity = 1;
	newArc->next = NULL;
	return newArc;
}

void output_arcGVZ ( char * outfile, boolean IsContig )
{
	ARC * pArc;
	preARC * pPreArc;
	char name[256];
	FILE * fp;
	unsigned int i;
	sprintf ( name, "%s.arc.gvz", outfile );
	fp = ckopen ( name, "w" );
	fprintf ( fp, "digraph G{\n" );
	fprintf ( fp, "\tsize=\"512,512\";\n" );

	for ( i = 1; i <= num_ed; i++ )
	{
		if ( IsContig )
		{
			pPreArc = contig_array[i].arcs;

			while ( pPreArc )
			{
				fprintf ( fp, "\tC%d -> C%d[label =\"%d\"];\n", i, pPreArc->to_ed, pPreArc->multiplicity );
				pPreArc = pPreArc->next;
			}
		}
		else
		{
			pArc = edge_array[i].arcs;

			while ( pArc )
			{
				fprintf ( fp, "\tC%d -> C%d[label =\"%d\"];\n", i, pArc->to_ed, pArc->multiplicity );
				pArc = pArc->next;
			}
		}
	}

	fprintf ( fp, "}\n" );
	fclose ( fp );
}

/**************** below this line all codes are about ARC ****************/
#define ARCBLOCKSIZE 100000
void createArcMemo ()
{
	if ( !arc_mem_manager )
	{
		arc_mem_manager = createMem_manager ( ARCBLOCKSIZE, sizeof ( ARC ) );
	}
	else
	{
		fprintf ( stderr, "Warning from createArcMemo: arc_mem_manager is a active pointer.\n" );
	}
}

void destroyArcMem ()
{
	freeMem_manager ( arc_mem_manager );
	arc_mem_manager = NULL;
}

ARC * allocateArc ( unsigned int edgeid )
{
	arcCounter++;
	ARC * newArc;
	newArc = ( ARC * ) getItem ( arc_mem_manager );
	newArc->to_ed = edgeid;
	newArc->multiplicity = 1;
	newArc->prev = NULL;
	newArc->next = NULL;
	return newArc;
}

void dismissArc ( ARC * arc )
{
	returnItem ( arc_mem_manager, arc );
}

/***************** below this line all codes are about lookup table *****************/

void createArcLookupTable ()
{
	if ( !arcLookupTable )
	{
		arcLookupTable = ( ARC ** ) ckalloc ( ( 3 * num_ed + 1 ) * sizeof ( ARC * ) );
	}
}

void deleteArcLookupTable ()
{
	if ( arcLookupTable )
	{
		free ( ( void * ) arcLookupTable );
		arcLookupTable = NULL;
	}
}

void putArc2LookupTable ( unsigned int from_ed, ARC * arc )
{
	if ( !arc || !arcLookupTable )
	{
		return;
	}

	unsigned int index = 2 * from_ed + arc->to_ed;
	arc->nextInLookupTable = arcLookupTable[index];
	arcLookupTable[index] = arc;
}

static ARC * getArcInLookupTable ( unsigned int from_ed, unsigned int to_ed )
{
	unsigned int index = 2 * from_ed + to_ed;
	ARC * ite_arc = arcLookupTable[index];

	while ( ite_arc )
	{
		if ( ite_arc->to_ed == to_ed )
		{
			return ite_arc;
		}

		ite_arc = ite_arc->nextInLookupTable;
	}

	return NULL;
}

void removeArcInLookupTable ( unsigned int from_ed, unsigned int to_ed )
{
	unsigned int index = 2 * from_ed + to_ed;
	ARC * ite_arc = arcLookupTable[index];
	ARC * arc;

	if ( !ite_arc )
	{
		fprintf ( stderr, "RemoveArcInLookupTable: not found A.\n" );
		return;
	}

	if ( ite_arc->to_ed == to_ed )
	{
		arcLookupTable[index] = ite_arc->nextInLookupTable;
		return;
	}

	while ( ite_arc->nextInLookupTable && ite_arc->nextInLookupTable->to_ed != to_ed )
	{
		ite_arc = ite_arc->nextInLookupTable;
	}

	if ( ite_arc->nextInLookupTable )
	{
		arc = ite_arc->nextInLookupTable;
		ite_arc->nextInLookupTable = arc->nextInLookupTable;
		return;
	}

	fprintf ( stderr, "RemoveArcInLookupTable: not found B.\n" );
	return;
}

void recordArcsInLookupTable ()
{
	unsigned int i;
	ARC * ite_arc;

	for ( i = 1; i <= num_ed; i++ )
	{
		ite_arc = edge_array[i].arcs;

		while ( ite_arc )
		{
			putArc2LookupTable ( i, ite_arc );
			ite_arc = ite_arc->next;
		}
	}
}

/*************************************************
Function:
    getArcBetween
Description:
    Gets the arc between the edge "from_ed" and the edge "to_ed", then return it.
Input:
    1. from_ed:     the origin edge
    2. to_ed:       the destination edge
Output:
    None.
Return:
    The arc between the edge "from_ed" and the edge "to_ed".
*************************************************/
ARC * getArcBetween ( unsigned int from_ed, unsigned int to_ed )
{
	ARC * parc;

	if ( arcLookupTable )
	{
		parc = getArcInLookupTable ( from_ed, to_ed );
		return parc;
	}

	parc = edge_array[from_ed].arcs;

	while ( parc )
	{
		if ( parc->to_ed == to_ed )
		{
			return parc;
		}

		parc = parc->next;
	}

	return parc;
}

/*************************************************
Function:
    validArcCount
Description:
    Compute the count of the valid arc, whose weight is larger than "cutoff", then return it.
Input:
    1. arc:         the head of the "preARC" linked list
    2. cutoff:      the minimum requirement for the weight of the valid arc
Output:
    None.
Return:
    The number of the valid arc.
*************************************************/
int validArcCount ( preARC * arc, int cutoff )
{
	int arc_num = 0;

	while ( arc )
	{
		if ( arc->multiplicity >= cutoff )
			{ ++arc_num; }

		arc = arc->next;
	}

	return arc_num;
}

/*************************************************
Function:
    maxArcWeight
Description:
    Gets the maximum weight in all arcs, then return it.
Input:
    1. arc:     the head of the "preARC" linked list
Output:
    None.
Return:
    The maximum weight.
*************************************************/
unsigned int maxArcWeight ( preARC * arc )
{
	unsigned int max = 0;

	while ( arc )
	{
		if ( arc->multiplicity > max )
			{ max = arc->multiplicity; }

		arc = arc->next;
	}

	return max;
}


/*
 * cutTip_graph.c
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

static int caseA, caseB, caseC, caseD, caseE;

/*************************************************
Function:
    destroyEdge
Description:
    Removes the edge and its arc.
Input:
    1. edgeid:  the edge index.
Output:
    None.
Return:
    None.
*************************************************/
void destroyEdge ( unsigned int edgeid )
{
	unsigned int bal_ed = getTwinEdge ( edgeid );
	ARC * arc;

	if ( bal_ed == edgeid )
	{
		edge_array[edgeid].length = 0;
		return;
	}

	arc = edge_array[edgeid].arcs;

	while ( arc )
	{
		arc->bal_arc->to_ed = 0;
		arc = arc->next;
	}

	arc = edge_array[bal_ed].arcs;

	while ( arc )
	{
		arc->bal_arc->to_ed = 0;
		arc = arc->next;
	}

	edge_array[edgeid].arcs = NULL;
	edge_array[bal_ed].arcs = NULL;
	edge_array[edgeid].length = 0;
	edge_array[bal_ed].length = 0;
	edge_array[edgeid].deleted = 1;
	edge_array[bal_ed].deleted = 1;
	//printf("Destroyed %d and %d\n",edgeid,bal_ed);
}

/*************************************************
Function:
    arcCount
Description:
    Computes the count of arc of the edge.
Input:
    1. edgeid:      the edge index
    2. num:         used to record the count of arc of the edge
Output:
    None.
Return:
    The first arc of the edge.
*************************************************/
ARC * arcCount ( unsigned int edgeid, unsigned int * num )
{
	ARC * arc;
	ARC * firstValidArc = NULL;
	unsigned int count = 0;
	arc = edge_array[edgeid].arcs;

	while ( arc )
	{
		if ( arc->to_ed > 0 )
		{
			count++;

			if ( count == 1 )
			{
				firstValidArc = arc;
			}
			else if ( count > 1 )
			{
				*num = count;
				return firstValidArc;
			}
		}

		arc = arc->next;
	}

	*num = count;
	return firstValidArc;
}

/*            multiplicity < multiCutoff
              ====  - ====  - ====
              length < lenCutoff
*/
/*************************************************
Function:
    removeWeakEdges
Description:
    Checks if the edge is short and its in-degree(or out-degree) is weak.
Input:
    1. lenCutoff:       the cutoff for the length of edge
    2. multiCutoff:     the cutoff for the weight of preArc
Output:
    None.
Return:
    None.
*************************************************/
void removeWeakEdges ( int lenCutoff, unsigned int multiCutoff )
{
	unsigned int bal_ed;
	unsigned int arcRight_n, arcLeft_n;
	ARC * arcLeft, *arcRight;
	unsigned int i;
	int counter = 0;
	int round = 1;
	fprintf ( stderr, "Start to destroy weak inner edges.\n" );
	counter = 1;

	while ( counter )
	{
		counter = 0;

		for ( i = 1; i <= num_ed; i++ )
		{
			if ( edge_array[i].deleted || edge_array[i].length == 0 || edge_array[i].length > lenCutoff || EdSameAsTwin ( i ) )
			{
				continue;
			}

			bal_ed = getTwinEdge ( i );
			arcRight = arcCount ( i, &arcRight_n );

			if ( arcRight_n > 1 || !arcRight || arcRight->multiplicity > multiCutoff )
			{
				continue;
			}

			arcLeft = arcCount ( bal_ed, &arcLeft_n );

			if ( arcLeft_n > 1 || !arcLeft || arcLeft->multiplicity > multiCutoff )
			{
				continue;
			}

			destroyEdge ( i );
			counter++;
		}

		fprintf ( stderr, "%d weak inner edge(s) destroyed in cycle %d.\n", counter, round++ );
	}

	removeDeadArcs ();
	/*
	   linearConcatenate();
	   compactEdgeArray();
	 */
}

/*
                cvg < covCutoff

        1.      ====  - ====  - ====

                                      ====
        2.      ====  - ====  <
                                      ====
                ====
        3.               > ====  - ====
                ====

                ====                 ====
        4.               > ====  <
                ====                 ====

                length < lenCutoff
*/
/*************************************************
Function:
    removeLowCovEdges
Description:
    Checks if the edge is short and its coverage is low.
Input:
    1. lenCutoff:       the cutoff for the length of edge.
    2. covCutoff:       the cutoff for the coverage of edge.
Output:
    None.
Return:
    None.
*************************************************/
void removeLowCovEdges ( int lenCutoff, unsigned short covCutoff, boolean last )
{
	unsigned int bal_ed;
	unsigned int arcRight_n, arcLeft_n;
	ARC * arcLeft, *arcRight;
	unsigned int i;
	int counter = 0;

	for ( i = 1; i <= num_ed; i++ )
	{
		if ( edge_array[i].deleted || edge_array[i].cvg == 0 || edge_array[i].cvg > covCutoff * 10 || edge_array[i].length >= lenCutoff || EdSameAsTwin ( i ) || edge_array[i].length == 0 )
		{
			continue;
		}

		bal_ed = getTwinEdge ( i );
		arcRight = arcCount ( i, &arcRight_n );
		arcLeft = arcCount ( bal_ed, &arcLeft_n );

		if ( arcLeft_n < 1 || arcRight_n < 1 )
		{
			continue;
		}

		destroyEdge ( i );
		counter++;
	}

	fprintf ( stderr, "%d inner edge(s) with coverage lower than or equal to %d destroyed.\n", counter, covCutoff );
	removeDeadArcs ();
	linearConcatenate ( 0, last );
	compactEdgeArray ();
}

/*************************************************
Function:
    isUnreliableTip
Description:
    Checks if the tips is suitable to cut.
Input:
    1. edgeid:      the edge index
    2. cutLen:      the length requirement for tips
    3. strict:          the choice for the pattern of cutting tips
Output:
    None.
Return:
    0 if the tip couldn't be cut.
*************************************************/
boolean isUnreliableTip ( unsigned int edgeid, int cutLen, boolean strict )
{
	unsigned int arcRight_n, arcLeft_n;
	unsigned int bal_ed;
	unsigned int currentEd = edgeid;
	int length = 0;
	unsigned int mult = 0;
	ARC * arc, *activeArc = NULL, *tempArc;

	if ( edgeid == 0 )
	{
		return 0;
	}

	bal_ed = getTwinEdge ( edgeid );

	if ( bal_ed == edgeid )
	{
		return 0;
	}

	arcCount ( bal_ed, &arcLeft_n );

	if ( arcLeft_n > 0 )
	{
		return 0;
	}

	while ( currentEd )
	{
		arcCount ( bal_ed, &arcLeft_n );
		tempArc = arcCount ( currentEd, &arcRight_n );

		if ( arcLeft_n > 1 || arcRight_n > 1 )
		{
			break;
		}

		length += edge_array[currentEd].length;

		if ( tempArc )
		{
			activeArc = tempArc;
			currentEd = activeArc->to_ed;
			bal_ed = getTwinEdge ( currentEd );
		}
		else
		{
			currentEd = 0;
		}
	}

	if ( length >= cutLen )
	{
		return 0;
	}

	if ( currentEd == 0 )
	{
		caseB++;
		return 1;
	}

	if ( !strict )
	{
		if ( arcLeft_n < 2 )
		{
			length += edge_array[currentEd].length;
		}

		if ( length >= cutLen )
		{
			return 0;
		}
		else
		{
			caseC++;
			return 1;
		}
	}

	if ( arcLeft_n < 2 )
	{
		return 0;
	}

	if ( !activeArc )
	{
		fprintf ( stderr, "No activeArc while checking edge %d.\n", edgeid );
	}

	if ( activeArc->multiplicity == 1 )
	{
		caseD++;
		return 1;
	}

	for ( arc = edge_array[bal_ed].arcs; arc != NULL; arc = arc->next )
		if ( arc->multiplicity > mult )
		{
			mult = arc->multiplicity;
		}

	if ( mult > activeArc->multiplicity )
	{
		caseE++;
	}

	return mult > activeArc->multiplicity;
}

boolean isUnreliableTip_strict ( unsigned int edgeid, int cutLen )
{
	unsigned int arcRight_n, arcLeft_n;
	unsigned int bal_ed;
	unsigned int currentEd = edgeid;
	int length = 0;
	unsigned int mult = 0;
	ARC * arc, *activeArc = NULL, *tempArc;

	if ( edgeid == 0 )
	{
		return 0;
	}

	bal_ed = getTwinEdge ( edgeid );

	if ( bal_ed == edgeid )
	{
		return 0;
	}

	arcCount ( bal_ed, &arcLeft_n );

	if ( arcLeft_n > 0 )
	{
		return 0;
	}

	while ( currentEd )
	{
		arcCount ( bal_ed, &arcLeft_n );
		tempArc = arcCount ( currentEd, &arcRight_n );

		if ( arcLeft_n > 1 || arcRight_n > 1 )
		{
			if ( arcLeft_n == 0 || length == 0 )
			{
				return 0;
			}
			else
			{
				break;
			}
		}

		length += edge_array[currentEd].length;

		if ( length >= cutLen )
		{
			return 0;
		}

		if ( tempArc )
		{
			activeArc = tempArc;
			currentEd = activeArc->to_ed;
			bal_ed = getTwinEdge ( currentEd );
		}
		else
		{
			currentEd = 0;
		}
	}

	if ( currentEd == 0 )
	{
		caseA++;
		return 1;
	}

	if ( !activeArc )
	{
		fprintf ( stderr, "No activeArc while checking edge %d.\n", edgeid );
	}

	if ( activeArc->multiplicity == 1 )
	{
		caseB++;
		return 1;
	}

	for ( arc = edge_array[bal_ed].arcs; arc != NULL; arc = arc->next )
		if ( arc->multiplicity > mult )
		{
			mult = arc->multiplicity;
		}

	if ( mult > activeArc->multiplicity )
	{
		caseC++;
	}

	return mult > activeArc->multiplicity;
}

/*************************************************
Function:
    removeDeadArcs
Description:
    Removes unvalid arcs.
Input:
    None.
Output:
    None.
Return:
    None.
*************************************************/
void removeDeadArcs ()
{
	unsigned int i, count = 0;
	ARC * arc, *arc_temp;

	for ( i = 1; i <= num_ed; i++ )
	{
		arc = edge_array[i].arcs;

		while ( arc )
		{
			arc_temp = arc;
			arc = arc->next;

			if ( arc_temp->to_ed == 0 )
			{
				count++;
				edge_array[i].arcs = deleteArc ( edge_array[i].arcs, arc_temp );
			}
		}
	}

	fprintf ( stderr, "%d dead arc(s) removed.\n", count );
}

/*************************************************
Function:
    cutTipsInGraph
Description:
    Cuts the short tips.
Input:
    1. cutLen:  the cutoff for the total length of the tips
    2. strict:      the pattern of the cutting tips
Output:
    None.
Return:
    None.
*************************************************/
void cutTipsInGraph ( int cutLen, boolean strict, boolean last )
{
	int flag = 1;
	unsigned int i;

	if ( !cutLen )
	{
		cutLen = 2 * overlaplen;
	}

	fprintf ( stderr, "\nStrict: %d, cutoff length: %d.\n", strict, cutLen );

	if ( strict )
	{
		linearConcatenate ( 0, last );
	}

	caseA = caseB = caseC = caseD = caseE = 0;
	int round = 1;

	while ( flag )
	{
		flag = 0;

		for ( i = 1; i <= num_ed; i++ )
		{
			if ( edge_array[i].deleted )
			{
				continue;
			}

			if ( isUnreliableTip ( i, cutLen, strict ) )
			{
				destroyEdge ( i );
				flag++;
			}
		}

		fprintf ( stderr, "%d tips cut in cycle %d.\n", flag, round++ );
	}

	removeDeadArcs ();

	if ( strict )
	{
		fprintf ( stderr, "Case A %d, B %d C %d D %d E %d.\n", caseA, caseB, caseC, caseD, caseE );
	}

	linearConcatenate ( 0, last );
	compactEdgeArray ();
}

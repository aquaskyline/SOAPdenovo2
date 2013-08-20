/*
 * splitReps.c
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

//Note: branches number is greater than 4 for sparse
static const int MAX_BRANCH_NUM = 1024; // the limit for the branch.
static unsigned int involved[2049]; //record the upstream and downstream edge
static unsigned int lefts[1024];        //record the upstream  edge
static unsigned int rights[1024];       //record the downstream edge
static unsigned char gothrough[1024][1024]; //record the path from upstream edges to downstream edges

/*************************************************
Function:
    interferingCheck
Description:
    1. Records the upsteam and downstream edges of current edge.
    2. Checks if any of the upstream/downstream edge (or its twin edge) is the upstream edge
       and the downstream edge at the same time.
Input:
    1. edgeno:      the current edge
    2. repTimes:        the in-degree (or out-degree)
Output:
    None.
Return:
    1 if the repeat is not suitable to solve.
*************************************************/
static boolean interferingCheck ( unsigned int edgeno, int repTimes )
{
	int i, j, t;
	unsigned int bal_ed;
	involved[0] = edgeno;
	i = 1;

	for ( j = 0; j < repTimes; j++ )
	{
		involved[i++] = lefts[j];
	}

	for ( j = 0; j < repTimes; j++ )
	{
		involved[i++] = rights[j];
	}

	for ( j = 0; j < i - 1; j++ )
		for ( t = j + 1; t < i; t++ )
			if ( involved[j] == involved[t] )
			{
				return 1;
			}

	for ( j = 0; j < i; j++ )
	{
		bal_ed = getTwinEdge ( involved[j] );

		for ( t = 0; t < i; t++ )
			if ( bal_ed == involved[t] )
			{
				return 1;
			}
	}

	return 0;
}

/*************************************************
Function:
    arcCounts
Description:
    Calculates the num of the downstream edges.
Input:
    1. edgeid:      the edge index
Output:
    1. num:     the number of downstream edges
Return:
    The first downstream arc.
*************************************************/
static ARC * arcCounts ( unsigned int edgeid, unsigned int * num )
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
		}

		if ( count == 1 )
		{
			firstValidArc = arc;
		}

		arc = arc->next;
	}

	*num = count;
	return firstValidArc;
}

/*************************************************
Function:
    readOnEdge
Description:
    Checks if the read ID exists in the marker array of a edge.
Input:
    1. readid:      the read ID
    2. edge:            the edge
Output:
    None.
Return:
    1 if the read ID exists in array.
*************************************************/
static boolean readOnEdge ( long long readid, unsigned int edge )
{
	int index;
	int markNum;
	long long * marklist;

	if ( edge_array[edge].markers )
	{
		markNum = edge_array[edge].multi;
		marklist = edge_array[edge].markers;
	}
	else
	{
		return 0;
	}

	for ( index = 0; index < markNum; index++ )
	{
		if ( readid == marklist[index] )
		{
			return 1;
		}
	}

	return 0;
}

/*************************************************
Function:
    cntByReads
Description:
    Checks if the read supports the path acrossing three edges.
Input:
    1. left:            the left edge
    2. middle:      the middle edge
    3. right:           the right edge
Output:
    None.
Return:
    1 if the read supports the path.
*************************************************/
static long long cntByReads ( unsigned int left, unsigned int middle, unsigned int right )
{
	int markNum;
	long long * marklist;

	if ( edge_array[left].markers )
	{
		markNum = edge_array[left].multi;
		marklist = edge_array[left].markers;
	}
	else
	{
		return 0;
	}

	int index;
	long long readid;

	/*
	   if(middle==8553)
	   printf("%d markers on %d\n",markNum,left);
	 */
	for ( index = 0; index < markNum; index++ )
	{
		readid = marklist[index];

		if ( readOnEdge ( readid, middle ) && readOnEdge ( readid, right ) )
		{
			return readid;
		}
	}

	return 0;
}

/*
        -       -
          > - <
        -       -
*/
/*************************************************
Function:
    solvable
Description:
    Checks if a edge is within a solvable repeat structure according to the following criterions:
    1. in-degree equals to the out-degree.
    2. the in-degree of the downstream edge is 1.
    3. the out-degree of the upsteam edge is 1.
    4. no interleaved read path between upstream and downstream edges.
Input:
    1. edgeno:      the edge
Output:
    None.
Return:
    0 if it is not.
*************************************************/
unsigned int solvable ( unsigned int edgeno )
{
	if ( EdSameAsTwin ( edgeno ) || edge_array[edgeno].multi == 255 )
	{
		return 0;
	}

	unsigned int bal_ed = getTwinEdge ( edgeno );
	unsigned int arcRight_n, arcLeft_n;
	unsigned int counter;
	unsigned int i, j;
	unsigned int branch, bal_branch;
	ARC * parcL, *parcR;
	parcL = arcCounts ( bal_ed, &arcLeft_n );

	if ( arcLeft_n < 2 )
	{
		return 0;
	}

	parcR = arcCounts ( edgeno, &arcRight_n );

	if ( arcLeft_n != arcRight_n )
	{
		return 0;
	}

	// check each right branch only has one upsteam connection
	/*
	   if(edgeno==2551){
	   for(i=0;i<arcLeft_n;i++)
	   printf("%d,",lefts[i]);
	   printf("__left to %d\n",edgeno);
	   for(j=0;j<arcRight_n;j++)
	   printf("%d,",rights[j]);
	   printf("__right to %d\n",edgeno);
	   }
	 */
	arcRight_n = 0;

	while ( parcR )
	{
		if ( parcR->to_ed == 0 )
		{
			parcR = parcR->next;
			continue;
		}

		branch = parcR->to_ed;

		if ( EdSameAsTwin ( branch ) || edge_array[branch].multi == 255 )
		{
			return 0;
		}

		rights[arcRight_n++] = branch;

		if ( arcRight_n >= MAX_BRANCH_NUM )
		{
			fprintf ( stderr, "ERROR: right arc number is bigger than the max %d.\n", MAX_BRANCH_NUM );
			exit ( -1 );
		}

		bal_branch = getTwinEdge ( branch );
		arcCounts ( bal_branch, &counter );

		if ( counter != 1 )
		{
			return 0;
		}

		parcR = parcR->next;
	}

	// check if each left branch only has one downsteam connection
	arcLeft_n = 0;

	while ( parcL )
	{
		if ( parcL->to_ed == 0 )
		{
			parcL = parcL->next;
			continue;
		}

		branch = parcL->to_ed;

		if ( EdSameAsTwin ( branch ) || edge_array[branch].multi == 255 )
		{
			return 0;
		}

		bal_branch = getTwinEdge ( branch );
		lefts[arcLeft_n++] = bal_branch;

		if ( arcLeft_n >= MAX_BRANCH_NUM )
		{
			fprintf ( stderr, "ERROR: left arc number is bigger than the max %d.\n", MAX_BRANCH_NUM );
			exit ( -1 );
		}

		arcCounts ( bal_branch, &counter );

		if ( counter != 1 )
		{
			return 0;
		}

		parcL = parcL->next;
	}

	//check if reads indicate one to one connection between upsteam and downstream edges

	for ( i = 0; i < arcLeft_n; i++ )
	{
		counter = 0;

		for ( j = 0; j < arcRight_n; j++ )
		{
			gothrough[i][j] = cntByReads ( lefts[i], edgeno, rights[j] ) == 0 ? 0 : 1;
			counter += gothrough[i][j];

			if ( counter > 1 )
			{
				return 0;
			}
		}

		if ( counter != 1 )
		{
			return 0;
		}
	}

	for ( j = 0; j < arcRight_n; j++ )
	{
		counter = 0;

		for ( i = 0; i < arcLeft_n; i++ )
		{
			counter += gothrough[i][j];
		}

		if ( counter != 1 )
		{
			return 0;
		}
	}

	return arcLeft_n;
}

/*************************************************
Function:
    cp1edge
Description:
    Copies the info of one edge to the other.
Input:
    1. source:  the source edge
    2. target:      the target edge
Output:
    None.
Return:
    The target edge.
*************************************************/
static unsigned int cp1edge ( unsigned int source, unsigned int target )
{
	int length = edge_array[source].length;
	char * tightSeq;
	int index;
	unsigned int bal_source = getTwinEdge ( source );
	unsigned int bal_target;

	if ( bal_source > source )
	{
		bal_target = target + 1;
	}
	else
	{
		bal_target = target;
		target = target + 1;
	}

	tightSeq = ( char * ) ckalloc ( ( length / 4 + 1 ) * sizeof ( char ) );

	for ( index = 0; index < length / 4 + 1; index++ )
	{
		tightSeq[index] = edge_array[source].seq[index];
	}

	edge_array[target].length = length;
	edge_array[target].cvg = edge_array[source].cvg;
	edge_array[target].to_vt = edge_array[source].to_vt;
	edge_array[target].from_vt = edge_array[source].from_vt;
	edge_array[target].seq = tightSeq;
	edge_array[target].bal_edge = edge_array[source].bal_edge;
	edge_array[target].rv = NULL;
	edge_array[target].arcs = NULL;
	edge_array[target].markers = NULL;
	edge_array[target].flag = 0;
	edge_array[target].deleted = 0;
	tightSeq = ( char * ) ckalloc ( ( length / 4 + 1 ) * sizeof ( char ) );

	for ( index = 0; index < length / 4 + 1; index++ )
	{
		tightSeq[index] = edge_array[bal_source].seq[index];
	}

	edge_array[bal_target].length = length;
	edge_array[bal_target].cvg = edge_array[bal_source].cvg;
	edge_array[bal_target].to_vt = edge_array[bal_source].to_vt;
	edge_array[bal_target].from_vt = edge_array[bal_source].from_vt;
	edge_array[bal_target].seq = tightSeq;
	edge_array[bal_target].bal_edge = edge_array[bal_source].bal_edge;
	edge_array[bal_target].rv = NULL;
	edge_array[bal_target].arcs = NULL;
	edge_array[bal_target].markers = NULL;
	edge_array[bal_target].flag = 0;
	edge_array[bal_target].deleted = 0;
	return target;
}

/*************************************************
Function:
    moveArc2cp
Description:
    Moves the upstream and downstream arc information of source edge to target edge.
Input:
    1. leftEd:          upstream edge
    2. rightEd:     dowmstream edge
    3. source:      the source edge
    4. target:          the target edge
Output:
    None.
Return:
    None.
*************************************************/
static void moveArc2cp ( unsigned int leftEd, unsigned int rightEd, unsigned int source, unsigned int target )
{
	unsigned int bal_left = getTwinEdge ( leftEd );
	unsigned int bal_right = getTwinEdge ( rightEd );
	unsigned int bal_source = getTwinEdge ( source );
	unsigned int bal_target = getTwinEdge ( target );
	ARC * arc;
	ARC * newArc, *twinArc;
	//between left and source
	arc = getArcBetween ( leftEd, source );
	arc->to_ed = 0;
	newArc = allocateArc ( target );
	newArc->multiplicity = arc->multiplicity;
	newArc->prev = NULL;
	newArc->next = edge_array[leftEd].arcs;

	if ( edge_array[leftEd].arcs )
	{
		edge_array[leftEd].arcs->prev = newArc;
	}

	edge_array[leftEd].arcs = newArc;
	arc = getArcBetween ( bal_source, bal_left );
	arc->to_ed = 0;
	twinArc = allocateArc ( bal_left );
	twinArc->multiplicity = arc->multiplicity;
	twinArc->prev = NULL;
	twinArc->next = NULL;
	edge_array[bal_target].arcs = twinArc;
	newArc->bal_arc = twinArc;
	twinArc->bal_arc = newArc;
	//between source and right
	arc = getArcBetween ( source, rightEd );
	arc->to_ed = 0;
	newArc = allocateArc ( rightEd );
	newArc->multiplicity = arc->multiplicity;
	newArc->prev = NULL;
	newArc->next = NULL;
	edge_array[target].arcs = newArc;
	arc = getArcBetween ( bal_right, bal_source );
	arc->to_ed = 0;
	twinArc = allocateArc ( bal_target );
	twinArc->multiplicity = arc->multiplicity;
	twinArc->prev = NULL;
	twinArc->next = edge_array[bal_right].arcs;

	if ( edge_array[bal_right].arcs )
	{
		edge_array[bal_right].arcs->prev = twinArc;
	}

	edge_array[bal_right].arcs = twinArc;
	newArc->bal_arc = twinArc;
	twinArc->bal_arc = newArc;
}

/*************************************************
Function:
    split1edge
Description:
    Copies the center edge to a new edge and updates related arc information if the center edge is in a solvable repeat sturcture.
Input:
    1. edgeno:      the edge index , which is the repeat
    2. repTimes:        the in-degree( or out-degree)
Output:
    None.
Return:
    None.
*************************************************/
static void split1edge ( unsigned int edgeno, int repTimes )
{
	int i, j;
	unsigned int target;

	for ( i = 1; i < repTimes; i++ )
	{
		for ( j = 0; j < repTimes; j++ )
			if ( gothrough[i][j] > 0 )  // a path supported by read
			{
				break;
			}

		target = cp1edge ( edgeno, extraEdgeNum );
		moveArc2cp ( lefts[i], rights[j], edgeno, target );
		extraEdgeNum += 2;
	}
}

static void debugging ( unsigned int i )
{
	ARC * parc;
	parc = edge_array[i].arcs;

	if ( !parc )
	{
		fprintf ( stderr, "No downward connection for %d.\n", i );
	}

	while ( parc )
	{
		fprintf ( stderr, "%d -> %d\n", i, parc->to_ed );
		parc = parc->next;
	}
}

/*************************************************
Function:
    solveReps
Description:
    Uses the read path information to solve the repeats.
Input:
    None.
Output:
    None.
Return:
    None.
*************************************************/
void solveReps ()
{
	unsigned int i;
	unsigned int repTime;
	int counter = 0;
	boolean flag;
	//debugging(30514);
	extraEdgeNum = num_ed + 1;

	for ( i = 1; i <= num_ed; i++ )
	{
		repTime = solvable ( i );

		if ( repTime == 0 )
		{
			continue;
		}

		flag = interferingCheck ( i, repTime );

		if ( flag )
		{
			continue;
		}

		split1edge ( i, repTime );
		counter++;  //+= 2*(repTime-1);

		if ( EdSmallerThanTwin ( i ) )
		{
			i++;
		}
	}

	fprintf ( stderr, "%d repeat(s) are solvable, %d more edge(s).\n", counter, extraEdgeNum - 1 - num_ed );
	num_ed = extraEdgeNum - 1;
	removeDeadArcs ();

	if ( markersArray )
	{
		free ( ( void * ) markersArray );
		markersArray = NULL;
	}
}

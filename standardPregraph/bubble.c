/*
 * bubble.c
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
#include "dfibHeap.h"
#include "fibHeap.h"

#define false 0
#define true 1

#define SLOW_TO_FAST 1
#define FAST_TO_SLOW 0

#define MAXREADLENGTH 100
#define MAXCONNECTION 100

static int MAXNODELENGTH;   // the limit for the edge in the path
static int DIFF;        // the mininum for  the difference between the paths

static unsigned int outNodeArray[MAXCONNECTION];    //
static ARC * outArcArray[MAXCONNECTION];
static boolean HasChanged;  // whether reset the arc

static const int INDEL = 0;
static const int SIM[4][4] =    // the score matrix of comparison
{
	{1, 0, 0, 0},
	{0, 1, 0, 0},
	{0, 0, 1, 0},
	{0, 0, 0, 1}
};

//static variables
static READINTERVAL * fastPath;     // used to record the ordered edges, which is the saved path
static READINTERVAL * slowPath; // used to record the ordered edges, which is the merged path

static char fastSequence[MAXREADLENGTH];        // used to record the sequence of the fast path
static char slowSequence[MAXREADLENGTH];        // used to record the sequence of the slow path

static int fastSeqLength;       // the length of the sequence of the fast path
static int slowSeqLength;       // the length of the sequence of the slow path

static Time * times; // record the weight from the upstream edge to the current edge and  used to decide which upstream edge is better
static unsigned int * previous; // record the upstream edge
static unsigned int expCounter;
static unsigned int * expanded;
static double cutoff;   // the mini difference between the paths

static int Fmatrix[MAXREADLENGTH + 1][MAXREADLENGTH + 1];   //the score matrix of comparing the paths
static int slowToFastMapping[MAXREADLENGTH + 1];        // the edge in the slow path map to the fast path
static int fastToSlowMapping[MAXREADLENGTH + 1];        // the edge in the fast path map to the slow path

static DFibHeapNode ** dheapNodes;
static DFibHeap * dheap;

static unsigned int activeNode;

//static ARC *activeArc;
static unsigned int startingNode;
static int progress;

static unsigned int * eligibleStartingPoints;

// DEBUG
static long long caseA, caseB, caseC, caseD, caseE;
static long long dnodeCounter;
static long long rnodeCounter;
static long long btCounter;
static long long cmpCounter;
static long long simiCounter;
//static long long pinCounter;
static long long replaceCounter;
static long long getArcCounter;

// END OF DEBUG

/*
static void output_contig1(int id, EDGE *edge)
{
    int i;
    Kmer kmer;
    char ch;
    char kmerSeq[100];
    printf(">%d_\n",id);
    kmer = vt_array[edge->from_vt].kmer;
    for(i=overlaplen-1;i>=0;i--){
        ch = kmer&3;
        kmer >>= 2;
        kmerSeq[i] = ch;
    }
    for(i=0;i<overlaplen;i++)
        printf("%c",int2base((int)kmerSeq[i]));

    for(i=0;i<edge->length;i++){
    printf("%c",int2base((int)getCharInTightString(edge->seq,i)));
    if((i+overlaplen+1)%100==0)
            printf("\n");
    }
    printf("\n");
}*/
static void output_seq ( char * seq, int length, FILE * fp, unsigned int from_vt, unsigned int dest )
{
	int i;
	Kmer kmer;
	kmer = vt_array[from_vt].kmer;
	printKmerSeq ( fp, kmer );
	fprintf ( fp, " " );

	for ( i = 0; i < length; i++ )
	{
		fprintf ( fp, "%c", int2base ( ( int ) seq[i] ) );
	}

	if ( edge_array[dest].seq )
	{
		fprintf ( fp, " %c\n", int2base ( ( int ) getCharInTightString ( edge_array[dest].seq, 0 ) ) );
	}
	else
	{
		fprintf ( fp, " N\n" );
	}
}

static void print_path ( FILE * fp )
{
	READINTERVAL * marker;
	marker = fastPath->nextInRead;

	while ( marker->nextInRead )
	{
		fprintf ( fp, "%u ", marker->edgeid );
		marker = marker->nextInRead;
	}

	fprintf ( fp, "\n" );
	marker = slowPath->nextInRead;

	while ( marker->nextInRead )
	{
		fprintf ( fp, "%u ", marker->edgeid );
		marker = marker->nextInRead;
	}

	fprintf ( fp, "\n" );
}

static void output_pair ( int lengthF, int lengthS, FILE * fp, int nodeF, int nodeS, boolean merged, unsigned int source, unsigned int destination )
{
	fprintf ( fp, "$$ %d vs %d $$ %d\n", nodeF, nodeS, merged );
	output_seq ( fastSequence, lengthF, fp, edge_array[source].to_vt, destination );
	output_seq ( slowSequence, lengthS, fp, edge_array[source].to_vt, destination );
}

/*************************************************
Function:
    resetNodeStatus
Description:
    Resets the status of the edge.
Input:
    None.
Output:
    None.
Return:
    None.
*************************************************/
static void resetNodeStatus ()
{
	unsigned int index;
	ARC * arc;
	unsigned int bal_ed;

	for ( index = 1; index <= num_ed; index++ )
	{
		if ( EdSameAsTwin ( index ) )
		{
			edge_array[index].multi = 1;
			continue;
		}

		arc = edge_array[index].arcs;
		bal_ed = getTwinEdge ( index );

		while ( arc )
		{
			if ( arc->to_ed == bal_ed )
			{
				break;
			}

			arc = arc->next;
		}

		if ( arc )
		{
			edge_array[index].multi = 1;
			edge_array[bal_ed].multi = 1;
			index++;
			continue;
		}

		arc = edge_array[bal_ed].arcs;

		while ( arc )
		{
			if ( arc->to_ed == index )
			{
				break;
			}

			arc = arc->next;
		}

		if ( arc )
		{
			edge_array[index].multi = 1;
			edge_array[bal_ed].multi = 1;
		}
		else
		{
			edge_array[index].multi = 0;
			edge_array[bal_ed].multi = 0;
		}

		index++;
	}
}

/*
static void determineEligibleStartingPoints()
{
    long long index,counter=0;
    unsigned int node;
    unsigned int maxmult;
    ARC *parc;
    FibHeap *heap = newFibHeap();

    for(index=1;index<=num_ed;index++){
        if(edge_array[index].deleted||edge_array[index].length<1)
            continue;
        maxmult = counter = 0;
        parc = edge_array[index].arcs;
        while(parc){
            if(parc->multiplicity > maxmult)
                maxmult = parc->multiplicity;

            parc = parc->next;
        }
        if(maxmult<1){
            continue;
        }
        insertNodeIntoHeap(heap,-maxmult,index);
    }
    counter = 0;
    while((index=removeNextNodeFromHeap(heap))!=0){
        eligibleStartingPoints[counter++] = index;
    }

    destroyHeap(heap);
    printf("%lld edges out of %d are eligible starting points\n",counter,num_ed);
}
*/

/*************************************************
Function:
    nextStartingPoint
Description:
    Gets the next start edge.
Input:
    None.
Output:
    None.
Return:
    The next start edge.
*************************************************/
static unsigned int nextStartingPoint ()
{
	unsigned int index = 1;
	unsigned int result = 0;

	for ( index = progress + 1; index < num_ed; index++ )
	{
		result = index;

		if ( edge_array[index].deleted || edge_array[index].length < 1 )
		{
			continue;
		}

		if ( result == 0 )
		{
			return 0;
		}

		if ( edge_array[result].multi > 0 )
		{
			continue;
		}

		progress = index;
		return result;
	}

	return 0;
}

/*************************************************
Function:
    updateNodeStatus
Description:
    Updates the status of the edge.
Input:
    None.
Output:
    None.
Return:
    None.
*************************************************/
static void updateNodeStatus ()
{
	unsigned int i, node;

	for ( i = 0; i < expCounter; i++ )
	{
		node = expanded[i];
		edge_array[node].multi = 1;
		edge_array[getTwinEdge ( node )].multi = 1;
	}
}

/*************************************************
Function:
    getNodePrevious
Description:
    Gets the previous edge of the current edge , then return it.
Input:
    1. node:        the current edge.
Output:
    None.
Return:
    The previous edge of the current edge.
*************************************************/
unsigned int getNodePrevious ( unsigned int node )
{
	return previous[node];
}

/*************************************************
Function:
    isPreviousToNode
Description:
    Checks if the edge "previous" is the previous of the edge "target", if it is, return 1.
Input:
    1. previous:            the edge index
    2. target:          the edge index
Output:
    None.
Return:
    1 if the edge "previous" is the previous of the edge "target".
*************************************************/
static boolean isPreviousToNode ( unsigned int previous, unsigned int target )
{
	unsigned int currentNode = target;
	unsigned int previousNode = 0;
	Time targetTime = times[target];

	while ( currentNode )
	{
		if ( currentNode == previous )
		{
			return 1;
		}

		if ( currentNode == previousNode )
		{
			return 0;
		}

		if ( times[currentNode] != targetTime )
		{
			return 0;
		}

		previousNode = currentNode;
		currentNode = getNodePrevious ( currentNode );
	}

	return 0;
}

/*************************************************
Function:
    copySeq
Description:
    Copies the sequence from "sourseS" to "targetS".
Input:
    1. targetS:     the sequence
    2. sourceS:     the sequence
    3. pos:         the start position of the targetS
    4. length:      the length of the sourceS
Output:
    None.
Return:
    None.
*************************************************/
static void copySeq ( char * targetS, char * sourceS, int pos, int length )
{
	char ch;
	int i, index;
	index = pos;

	for ( i = 0; i < length; i++ )
	{
		ch = getCharInTightString ( sourceS, i );
		targetS[index++] = ch;
	}
}

/*************************************************
Function:
    extractSequence
Description:
    Copies all the sequence of the path to sequence.
Input:
    1. path:        a path consists of the ordered edge
    2. sequence:        used to record the sequence of the path
Output:
    None.
Return:
    The length of sequence.
*************************************************/
static int extractSequence ( READINTERVAL * path, char * sequence )
{
	READINTERVAL * marker;
	int seqLength, writeIndex;
	seqLength = writeIndex = 0;
	path->start = -10;
	marker = path->nextInRead;

	while ( marker->nextInRead )
	{
		marker->start = seqLength;
		seqLength += edge_array[marker->edgeid].length;
		marker = marker->nextInRead;
	}

	marker->start = seqLength;

	if ( seqLength > MAXREADLENGTH )
	{
		return 0;
	}

	marker = path->nextInRead;

	while ( marker->nextInRead )
	{
		if ( edge_array[marker->edgeid].length && edge_array[marker->edgeid].seq )
		{
			copySeq ( sequence, edge_array[marker->edgeid].seq, writeIndex, edge_array[marker->edgeid].length );
			writeIndex += edge_array[marker->edgeid].length;
		}

		/*
		   else if(edge_array[marker->edgeid].length==0)
		   printf("node %d with length 0 in this path\n",marker->edgeid);
		   else if(edge_array[marker->edgeid].seq==NULL)
		   printf("node %d without seq in this path\n",marker->edgeid);
		 */
		marker = marker->nextInRead;
	}

	return seqLength;
}

static int max ( int A, int B, int C )
{
	A = A >= B ? A : B;
	return ( A >= C ? A : C );
}

/*************************************************
Function:
    compareSequences
Description:
    Checks if the sequences are long enough and high similar. if not, return 0.
Input:
    1. sequence1:       the first sequence
    2. sequence2:       the second sequence
    3. length1:     the length of the first sequence
    4. length2:     the length of the second sequence
Output:
    None.
Return:
    0 if the bubble is not suitable to merge.
*************************************************/
static boolean compareSequences ( char * sequence1, char * sequence2, int length1, int length2 )
{
	int i, j;
	int maxLength;
	int Choice1, Choice2, Choice3;
	int maxScore;

	if ( length1 == 0 || length2 == 0 )
	{
		caseA++;
		return 0;
	}

	if ( abs ( ( int ) length1 - ( int ) length2 ) > 2 )
	{
		caseB++;
		return 0;
	}

	if ( length1 < overlaplen - 1 || length2 < overlaplen - 1 )
	{
		caseE++;
		return 0;
	}

	/*
	   if (length1 < overlaplen || length2 < overlaplen){
	   if(abs((int)length1 - (int)length2) > 3){
	   caseB++;
	   return 0;
	   }
	   }
	 */
	for ( i = 0; i <= length1; i++ )
	{
		Fmatrix[i][0] = 0;
	}

	for ( j = 0; j <= length2; j++ )
	{
		Fmatrix[0][j] = 0;
	}

	for ( i = 1; i <= length1; i++ )
	{
		for ( j = 1; j <= length2; j++ )
		{
			Choice1 = Fmatrix[i - 1][j - 1] + SIM[ ( int ) sequence1[i - 1]][ ( int ) sequence2[j - 1]];
			Choice2 = Fmatrix[i - 1][j] + INDEL;
			Choice3 = Fmatrix[i][j - 1] + INDEL;
			Fmatrix[i][j] = max ( Choice1, Choice2, Choice3 );
		}
	}

	maxScore = Fmatrix[length1][length2];
	maxLength = ( length1 > length2 ? length1 : length2 );

	if ( maxScore < maxLength - DIFF )
	{
		caseC++;
		return 0;
	}

	if ( ( 1 - ( double ) maxScore / maxLength ) > cutoff )
	{
		caseD++;
		return 0;
	}

	return 1;
}

/*************************************************
Function:
    mapSlowOntoFast
Description:
    Maps the sequence of the slow path to the one of the fast path.
Input:
    None.
Output:
    None.
Return:
    None.
*************************************************/
static void mapSlowOntoFast ()
{
	int slowIndex = slowSeqLength;
	int fastIndex = fastSeqLength;
	int fastn, slown;

	if ( slowIndex == 0 )
	{
		slowToFastMapping[0] = fastIndex;

		while ( fastIndex >= 0 )
		{
			fastToSlowMapping[fastIndex--] = 0;
		}

		return;
	}

	if ( fastIndex == 0 )
	{
		while ( slowIndex >= 0 )
		{
			slowToFastMapping[slowIndex--] = 0;
		}

		fastToSlowMapping[0] = slowIndex;
		return;
	}

	while ( slowIndex > 0 && fastIndex > 0 )
	{
		fastn = ( int ) fastSequence[fastIndex - 1]; //getCharInTightString(fastSequence,fastIndex-1);
		slown = ( int ) slowSequence[slowIndex - 1]; //getCharInTightString(slowSequence,slowIndex-1);

		if ( Fmatrix[fastIndex][slowIndex] == Fmatrix[fastIndex - 1][slowIndex - 1] + SIM[fastn][slown] )
		{
			fastToSlowMapping[--fastIndex] = --slowIndex;
			slowToFastMapping[slowIndex] = fastIndex;
		}
		else if ( Fmatrix[fastIndex][slowIndex] == Fmatrix[fastIndex - 1][slowIndex] + INDEL )
		{
			fastToSlowMapping[--fastIndex] = slowIndex - 1;
		}
		else if ( Fmatrix[fastIndex][slowIndex] == Fmatrix[fastIndex][slowIndex - 1] + INDEL )
		{
			slowToFastMapping[--slowIndex] = fastIndex - 1;
		}
		else
		{
			fprintf ( stderr, "Error in the step:  map the slow path to the fast path.\n" );
			abort ();
		}
	}

	while ( slowIndex > 0 )
	{
		slowToFastMapping[--slowIndex] = -1;
	}

	while ( fastIndex > 0 )
	{
		fastToSlowMapping[--fastIndex] = -1;
	}

	slowToFastMapping[slowSeqLength] = fastSeqLength;
	fastToSlowMapping[fastSeqLength] = slowSeqLength;
}

/*************************************************
Function:
    deleteArc
Description:
    Removes an arc from the double linked list and return the updated list
Input:
    1. arc_list:    the linked of the arc
    2. arc:     the deleted arc
Output:
    None.
Return:
    The new linked of the arc.
*************************************************/
ARC * deleteArc ( ARC * arc_list, ARC * arc )
{
	if ( arc->prev )
	{
		arc->prev->next = arc->next;
	}
	else
	{
		arc_list = arc->next;
	}

	if ( arc->next )
	{
		arc->next->prev = arc->prev;
	}

	/*
	   if(checkActiveArc&&arc==activeArc){
	   activeArc = arc->next;
	   }
	 */
	dismissArc ( arc );
	return arc_list;
}

/*************************************************
Function:
    addRv
Description:
    Add a rv to the head of a rv list.
Input:
    1. rv_list:     the linked of the path
    2. rv:          the new one, which will be added into the path
Output:
    None.
Return:
    The new linked of the path.
*************************************************/
static READINTERVAL * addRv ( READINTERVAL * rv_list, READINTERVAL * rv )
{
	rv->prevOnEdge = NULL;
	rv->nextOnEdge = rv_list;

	if ( rv_list )
	{
		rv_list->prevOnEdge = rv;
	}

	rv_list = rv;
	return rv_list;
}

/*************************************************
Function:
    deleteRv
Description:
    Removes a rv from the double linked list and return the updated list
Input:
    1. rv_list:     the linked of the path
    2. rv:          the deleted head from the path
Output:
    None.
Return:
    The new linked of the path.
*************************************************/
static READINTERVAL * deleteRv ( READINTERVAL * rv_list, READINTERVAL * rv )
{
	if ( rv->prevOnEdge )
	{
		rv->prevOnEdge->nextOnEdge = rv->nextOnEdge;
	}
	else
	{
		rv_list = rv->nextOnEdge;
	}

	if ( rv->nextOnEdge )
	{
		rv->nextOnEdge->prevOnEdge = rv->prevOnEdge;
	}

	return rv_list;
}

/*
static void disconnect(unsigned int from_ed, unsigned int to_ed)
{
    READINTERVAL *rv_temp;

    rv_temp = edge_array[from_ed].rv;
    while(rv_temp){
        if(!rv_temp->nextInRead||rv_temp->nextInRead->edgeid!=to_ed){
            rv_temp = rv_temp->nextOnEdge;
            continue;
        }
        rv_temp->nextInRead->prevInRead = NULL;
        rv_temp->nextInRead = NULL;

        rv_temp = rv_temp->nextOnEdge;
    }
}
*/

/*************************************************
Function:
    mapDistancesOntoPaths
Description:
    Sets the distance of the edge on the path.
Input:
    None.
Output:
    None.
Return:
    The total distance of the fast path.
*************************************************/
static int mapDistancesOntoPaths ()
{
	READINTERVAL * marker;
	int totalDistance = 0;
	marker = slowPath;

	while ( marker->nextInRead )
	{
		marker = marker->nextInRead;
		marker->start = totalDistance;
		totalDistance += edge_array[marker->edgeid].length;
		marker->bal_rv->start = totalDistance;
	}

	totalDistance = 0;
	marker = fastPath;

	while ( marker->nextInRead )
	{
		marker = marker->nextInRead;
		marker->start = totalDistance;
		totalDistance += edge_array[marker->edgeid].length;
		marker->bal_rv->start = totalDistance;
	}

	return totalDistance;
}


/*************************************************
Function:
    attachPath
Description:
    Attaches a path to the graph and mean while make the reverse complementary path of it
Input:
    1. path:        the target path
Output:
    None.
Return:
    None.
*************************************************/
static void attachPath ( READINTERVAL * path )
{
	READINTERVAL * marker, *bal_marker;
	unsigned int ed, bal_ed;
	marker = path;

	while ( marker )
	{
		ed = marker->edgeid;
		edge_array[ed].rv = addRv ( edge_array[ed].rv, marker );
		bal_ed = getTwinEdge ( ed );
		bal_marker = allocateRV ( -marker->readid, bal_ed );
		edge_array[bal_ed].rv = addRv ( edge_array[bal_ed].rv, bal_marker );

		if ( marker->prevInRead )
		{
			marker->prevInRead->bal_rv->prevInRead = bal_marker;
			bal_marker->nextInRead = marker->prevInRead->bal_rv;
		}

		bal_marker->bal_rv = marker;
		marker->bal_rv = bal_marker;
		marker = marker->nextInRead;
	}
}
static void detachPathSingle ( READINTERVAL * path )
{
	READINTERVAL * marker, *nextMarker;
	unsigned int ed;
	marker = path;

	while ( marker )
	{
		nextMarker = marker->nextInRead;
		ed = marker->edgeid;
		edge_array[ed].rv = deleteRv ( edge_array[ed].rv, marker );
		dismissRV ( marker );
		marker = nextMarker;
	}
}

/*************************************************
Function:
    detachPath
Description:
    Removes the path.
Input:
    1. path:        the linked of the path
Output:
    None.
Return:
    None.
*************************************************/
static void detachPath ( READINTERVAL * path )
{
	READINTERVAL * marker, *bal_marker, *nextMarker;
	unsigned int ed, bal_ed;
	marker = path;

	while ( marker )
	{
		nextMarker = marker->nextInRead;
		bal_marker = marker->bal_rv;
		ed = marker->edgeid;
		edge_array[ed].rv = deleteRv ( edge_array[ed].rv, marker );
		dismissRV ( marker );
		bal_ed = getTwinEdge ( ed );
		edge_array[bal_ed].rv = deleteRv ( edge_array[bal_ed].rv, bal_marker );
		dismissRV ( bal_marker );
		marker = nextMarker;
	}
}

static void remapNodeMarkersOntoNeighbour ( unsigned int source, unsigned int target )
{
	READINTERVAL * marker, *bal_marker;
	unsigned int bal_source = getTwinEdge ( source );
	unsigned int bal_target = getTwinEdge ( target );

	while ( ( marker = edge_array[source].rv ) != NULL )
	{
		edge_array[source].rv = deleteRv ( edge_array[source].rv, marker );
		marker->edgeid = target;
		edge_array[target].rv = addRv ( edge_array[target].rv, marker );
		bal_marker = marker->bal_rv;
		edge_array[bal_source].rv = deleteRv ( edge_array[bal_source].rv, bal_marker );
		bal_marker->edgeid = bal_target;
		edge_array[bal_target].rv = addRv ( edge_array[bal_target].rv, bal_marker );
	}
}

static void remapNodeInwardReferencesOntoNode ( unsigned int source, unsigned int target )
{
	ARC * arc;
	unsigned int destination;

	for ( arc = edge_array[source].arcs; arc != NULL; arc = arc->next )
	{
		destination = arc->to_ed;

		if ( destination == target || destination == source )
		{
			continue;
		}

		if ( previous[destination] == source )
		{
			previous[destination] = target;
		}
	}
}
static void remapNodeTimesOntoTargetNode ( unsigned int source, unsigned int target )
{
	Time nodeTime = times[source];
	unsigned int prevNode = previous[source];
	Time targetTime = times[target];

	if ( nodeTime == -1 )
	{
		return;
	}

	if ( prevNode == source )
	{
		times[target] = nodeTime;
		previous[target] = target;
	}
	else if ( targetTime == -1 || targetTime > nodeTime || ( targetTime == nodeTime && !isPreviousToNode ( target, prevNode ) ) )
	{
		times[target] = nodeTime;

		if ( prevNode != getTwinEdge ( source ) )
		{
			previous[target] = prevNode;
		}
		else
		{
			previous[target] = getTwinEdge ( target );
		}
	}

	remapNodeInwardReferencesOntoNode ( source, target );
	previous[source] = 0;
}

static void remapNodeTimesOntoNeighbour ( unsigned int source, unsigned int target )
{
	remapNodeTimesOntoTargetNode ( source, target );
	remapNodeTimesOntoTargetNode ( getTwinEdge ( source ), getTwinEdge ( target ) ); //questionable
}

static void destroyArc ( unsigned int from_ed, ARC * arc )
{
	unsigned int bal_dest;
	ARC * twinArc;

	if ( !arc )
	{
		return;
	}

	bal_dest = getTwinEdge ( arc->to_ed );
	twinArc = arc->bal_arc;
	removeArcInLookupTable ( from_ed, arc->to_ed );
	edge_array[from_ed].arcs = deleteArc ( edge_array[from_ed].arcs, arc );

	if ( bal_dest != from_ed )
	{
		removeArcInLookupTable ( bal_dest, getTwinEdge ( from_ed ) );
		edge_array[bal_dest].arcs = deleteArc ( edge_array[bal_dest].arcs, twinArc );
	}
}

static void createAnalogousArc ( unsigned int originNode, unsigned int destinationNode, ARC * refArc )
{
	ARC * arc, *twinArc;
	unsigned int destinationTwin;
	arc = getArcBetween ( originNode, destinationNode );

	if ( arc )
	{
		if ( refArc->bal_arc != refArc )
		{
			arc->multiplicity += refArc->multiplicity;
			arc->bal_arc->multiplicity += refArc->multiplicity;
		}
		else
		{
			arc->multiplicity += refArc->multiplicity / 2;
			arc->bal_arc->multiplicity += refArc->multiplicity / 2;
		}

		return;
	}

	arc = allocateArc ( destinationNode );
	arc->multiplicity = refArc->multiplicity;
	arc->prev = NULL;
	arc->next = edge_array[originNode].arcs;

	if ( edge_array[originNode].arcs )
	{
		edge_array[originNode].arcs->prev = arc;
	}

	edge_array[originNode].arcs = arc;
	putArc2LookupTable ( originNode, arc );
	destinationTwin = getTwinEdge ( destinationNode );

	if ( destinationTwin == originNode )
	{
		arc->bal_arc = arc;

		if ( refArc->bal_arc != refArc )
		{
			arc->multiplicity += refArc->multiplicity;
		}

		return;
	}

	twinArc = allocateArc ( getTwinEdge ( originNode ) );
	arc->bal_arc = twinArc;
	twinArc->bal_arc = arc;
	twinArc->multiplicity = refArc->multiplicity;
	twinArc->prev = NULL;
	twinArc->next = edge_array[destinationTwin].arcs;

	if ( edge_array[destinationTwin].arcs )
	{
		edge_array[destinationTwin].arcs->prev = twinArc;
	}

	edge_array[destinationTwin].arcs = twinArc;
	putArc2LookupTable ( destinationTwin, twinArc );
}

static void remapNodeArcsOntoTarget ( unsigned int source, unsigned int target )
{
	ARC * arc;

	if ( source == activeNode )
	{
		activeNode = target;
	}

	arc = edge_array[source].arcs;

	if ( !arc )
	{
		return;
	}

	while ( arc != NULL )
	{
		createAnalogousArc ( target, arc->to_ed, arc );
		destroyArc ( source, arc );
		arc = edge_array[source].arcs;
	}
}

static void remapNodeArcsOntoNeighbour ( unsigned int source, unsigned int target )
{
	remapNodeArcsOntoTarget ( source, target );
	remapNodeArcsOntoTarget ( getTwinEdge ( source ), getTwinEdge ( target ) );
}

static DFibHeapNode * getNodeDHeapNode ( unsigned int node )
{
	return dheapNodes[node];
}

static void setNodeDHeapNode ( unsigned int node, DFibHeapNode * dheapNode )
{
	dheapNodes[node] = dheapNode;
}

static void remapNodeFibHeapReferencesOntoNode ( unsigned int source, unsigned int target )
{
	DFibHeapNode * sourceDHeapNode = getNodeDHeapNode ( source );
	DFibHeapNode * targetDHeapNode = getNodeDHeapNode ( target );

	if ( sourceDHeapNode == NULL )
	{
		return;
	}

	if ( targetDHeapNode == NULL )
	{
		setNodeDHeapNode ( target, sourceDHeapNode );
		replaceValueInDHeap ( sourceDHeapNode, target );
	}
	else if ( getKey ( targetDHeapNode ) > getKey ( sourceDHeapNode ) )
	{
		setNodeDHeapNode ( target, sourceDHeapNode );
		replaceValueInDHeap ( sourceDHeapNode, target );
		destroyNodeInDHeap ( targetDHeapNode, dheap );
	}
	else
	{
		destroyNodeInDHeap ( sourceDHeapNode, dheap );
	}

	setNodeDHeapNode ( source, NULL );
}

static void combineCOV ( unsigned int source, int len_s, unsigned int target, int len_t )
{
	if ( len_s < 1 || len_t < 1 )
	{
		return;
	}

	int cov = ( len_s * edge_array[source].cvg + len_t * edge_array[target].cvg ) / len_t;
	edge_array[target].cvg = cov > MaxEdgeCov ? MaxEdgeCov : cov;
	edge_array[getTwinEdge ( target )].cvg = cov > MaxEdgeCov ? MaxEdgeCov : cov;
}

/*************************************************
Function:
    remapNodeOntoNeighbour
Description:
    Replaces the edge 'source' with the edge 'target' and adds the position info from 'source' to 'target'.
Input:
    1. source:      the edge index
    2. target:      the edge index
Output:
    None.
Return:
    None.
*************************************************/
static void remapNodeOntoNeighbour ( unsigned int source, unsigned int target )
{
	combineCOV ( source, edge_array[source].length, target, edge_array[target].length );
	remapNodeMarkersOntoNeighbour ( source, target );
	remapNodeTimesOntoNeighbour ( source, target ); //questionable
	remapNodeArcsOntoNeighbour ( source, target );
	remapNodeFibHeapReferencesOntoNode ( source, target );
	remapNodeFibHeapReferencesOntoNode ( getTwinEdge ( source ), getTwinEdge ( target ) );
	edge_array[source].deleted = 1;
	edge_array[getTwinEdge ( source )].deleted = 1;

	if ( startingNode == source )
	{
		startingNode = target;
	}

	if ( startingNode == getTwinEdge ( source ) )
	{
		startingNode = getTwinEdge ( target );
	}

	edge_array[source].length = 0;
	edge_array[getTwinEdge ( source )].length = 0;
}

static void connectInRead ( READINTERVAL * previous, READINTERVAL * next )
{
	if ( previous )
	{
		previous->nextInRead = next;

		if ( next )
		{
			previous->bal_rv->prevInRead = next->bal_rv;
		}
		else
		{
			previous->bal_rv->prevInRead = NULL;
		}
	}

	if ( next )
	{
		next->prevInRead = previous;

		if ( previous )
		{
			next->bal_rv->nextInRead = previous->bal_rv;
		}
		else
		{
			next->bal_rv->nextInRead = NULL;
		}
	}
}

static int remapBackOfNodeMarkersOntoNeighbour ( unsigned int source, READINTERVAL * sourceMarker, unsigned int target, READINTERVAL * targetMarker, boolean slowToFast )
{
	READINTERVAL * marker, *newMarker, *bal_new, *previousMarker;
	int halfwayPoint, halfwayPointOffset, breakpoint;
	int * targetToSourceMapping, *sourceToTargetMapping;
	unsigned int bal_ed;
	int targetFinish = targetMarker->bal_rv->start;
	int sourceStart = sourceMarker->start;
	int sourceFinish = sourceMarker->bal_rv->start;
	int alignedSourceLength = sourceFinish - sourceStart;
	int realSourceLength = edge_array[source].length;

	if ( slowToFast )
	{
		sourceToTargetMapping = slowToFastMapping;
		targetToSourceMapping = fastToSlowMapping;
	}
	else
	{
		sourceToTargetMapping = fastToSlowMapping;
		targetToSourceMapping = slowToFastMapping;
	}

	if ( alignedSourceLength > 0 && targetFinish > 0 )
	{
		halfwayPoint = targetToSourceMapping[targetFinish - 1] - sourceStart + 1;
		halfwayPoint *= realSourceLength;
		halfwayPoint /= alignedSourceLength;
	}
	else
	{
		halfwayPoint = 0;
	}

	if ( halfwayPoint < 0 )
	{
		halfwayPoint = 0;
	}

	if ( halfwayPoint > realSourceLength )
	{
		halfwayPoint = realSourceLength;
	}

	halfwayPointOffset = realSourceLength - halfwayPoint;
	bal_ed = getTwinEdge ( target );

	for ( marker = edge_array[source].rv; marker != NULL; marker = marker->nextOnEdge )
	{
		if ( marker->prevInRead && marker->prevInRead->edgeid == target )
		{
			continue;
		}

		newMarker = allocateRV ( marker->readid, target );
		edge_array[target].rv = addRv ( edge_array[target].rv, newMarker );
		bal_new = allocateRV ( -marker->readid, bal_ed );
		edge_array[bal_ed].rv = addRv ( edge_array[bal_ed].rv, bal_new );
		newMarker->bal_rv = bal_new;
		bal_new->bal_rv = newMarker;
		newMarker->start = marker->start;

		if ( realSourceLength > 0 )
		{
			breakpoint = halfwayPoint + marker->start;
		}
		else
		{
			breakpoint = marker->start;
		}

		bal_new->start = breakpoint;
		marker->start = breakpoint;
		previousMarker = marker->prevInRead;
		connectInRead ( previousMarker, newMarker );
		connectInRead ( newMarker, marker );
	}

	return halfwayPointOffset;
}

static void printKmer ( Kmer kmer )
{
	printKmerSeq ( stderr, kmer );
	fprintf ( stderr, "\n" );
}

static int splitNodeDescriptor ( unsigned int source, unsigned int target, int offset )
{
	int originalLength = edge_array[source].length;
	int backLength = originalLength - offset;
	int index, seqLen;
	char * tightSeq, nt, *newSeq;
	unsigned int bal_source = getTwinEdge ( source );
	unsigned int bal_target = getTwinEdge ( target );
	edge_array[source].length = offset;
	edge_array[bal_source].length = offset;
	edge_array[source].flag = 1;
	edge_array[bal_source].flag = 1;

	if ( target != 0 )
	{
		edge_array[target].length = backLength;
		edge_array[bal_target].length = backLength;
		free ( ( void * ) edge_array[target].seq );
		edge_array[target].seq = NULL;
		free ( ( void * ) edge_array[bal_target].seq );
		edge_array[bal_target].seq = NULL;
	}

	if ( backLength == 0 )
	{
		return 0;
	}

	tightSeq = edge_array[source].seq;
	seqLen = backLength / 4 + 1;

	if ( target != 0 )
	{
		edge_array[target].flag = 1;
		edge_array[bal_target].flag = 1;
		newSeq = ( char * ) ckalloc ( seqLen * sizeof ( char ) );
		edge_array[target].seq = newSeq;

		for ( index = 0; index < backLength; index++ )
		{
			nt = getCharInTightString ( tightSeq, index );
			writeChar2tightString ( nt, newSeq, index );
		}
	}

	//source node
	for ( index = backLength; index < originalLength; index++ )
	{
		nt = getCharInTightString ( tightSeq, index );
		writeChar2tightString ( nt, tightSeq, index - backLength );
	}

	if ( target == 0 )
	{
		return backLength;
	}

	//target twin
	tightSeq = edge_array[bal_source].seq;
	newSeq = ( char * ) ckalloc ( seqLen * sizeof ( char ) );
	edge_array[bal_target].seq = newSeq;

	for ( index = offset; index < originalLength; index++ )
	{
		nt = getCharInTightString ( tightSeq, index );
		writeChar2tightString ( nt, newSeq, index - offset );
	}

	return backLength;
}

static void remapBackOfNodeDescriptorOntoNeighbour ( unsigned int source, unsigned int target, boolean slowToFast, int offset )
{
	unsigned int bal_source = getTwinEdge ( source );
	unsigned int bal_target = getTwinEdge ( target );
	Kmer source_from_vt_kmer , bal_source_from_vt_kmer;
	Kmer word;
	int index;
	char nt;
	int backlength ;

	if ( slowToFast )
	{
		backlength = splitNodeDescriptor ( source, 0, offset );
		edge_array[source].from_vt = edge_array[target].to_vt;
		edge_array[bal_source].to_vt = edge_array[bal_target].from_vt;
	}
	else
	{
		backlength = splitNodeDescriptor ( source, target, offset );
		source_from_vt_kmer = vt_array[edge_array[source].from_vt].kmer;
		bal_source_from_vt_kmer = vt_array[edge_array[bal_source].to_vt].kmer;
		edge_array[target].from_vt = new_num_vt;

		if ( new_num_vt + 1 > num_kmer_limit )
		{
			fprintf ( stderr, "Error : Number of vertex is out of range.\n" );
			exit ( -1 );
		}

		vt_array[new_num_vt++].kmer = source_from_vt_kmer;
		edge_array[bal_target].to_vt = new_num_vt;

		if ( new_num_vt + 1 > num_kmer_limit )
		{
			fprintf ( stderr, "Error : Number of vertex is out of range.\n" );
			exit ( -1 );
		}

		vt_array[new_num_vt++].kmer = bal_source_from_vt_kmer;
		word = vt_array[edge_array[target].from_vt].kmer;

		for ( index = 0; index < backlength; index++ )
		{
			nt = getCharInTightString ( edge_array[target].seq, index );
			word = nextKmer ( word, nt );
		}

		edge_array[target].to_vt = new_num_vt;

		if ( new_num_vt + 1 > num_kmer_limit )
		{
			fprintf ( stderr, "Error : Number of vertex is out of range.\n" );
			exit ( -1 );
		}

		vt_array[new_num_vt++].kmer = word;
		edge_array[source].from_vt = new_num_vt;

		if ( new_num_vt + 1 > num_kmer_limit )
		{
			fprintf ( stderr, "Error : Number of vertex is out of range.\n" );
			exit ( -1 );
		}

		vt_array[new_num_vt++].kmer = word;
		word = vt_array[edge_array[bal_source].from_vt].kmer;

		for ( index = 0; index < offset; index++ )
		{
			nt = getCharInTightString ( edge_array[bal_source].seq, index );
			word = nextKmer ( word, nt );
		}

		edge_array[bal_target].from_vt = new_num_vt;

		if ( new_num_vt + 1 > num_kmer_limit )
		{
			fprintf ( stderr, "Error : Number of vertex is out of range.\n" );
			exit ( -1 );
		}

		vt_array[new_num_vt++].kmer = word;
		edge_array[bal_source].to_vt = new_num_vt;

		if ( new_num_vt + 1 > num_kmer_limit )
		{
			fprintf ( stderr, "Error : Number of vertex is out of range.\n" );
			exit ( -1 );
		}

		vt_array[new_num_vt++].kmer = word;
	}
}

static void remapBackOfNodeTimesOntoNeighbour ( unsigned int source, unsigned int target )
{
	Time targetTime = times[target];
	Time nodeTime = times[source];
	unsigned int twinTarget = getTwinEdge ( target );
	unsigned int twinSource = getTwinEdge ( source );
	unsigned int previousNode;

	if ( nodeTime != -1 )
	{
		previousNode = previous[source];

		if ( previousNode == source )
		{
			times[target] = nodeTime;
			previous[target] = target;
		}
		else if ( targetTime == -1 || targetTime > nodeTime || ( targetTime == nodeTime && !isPreviousToNode ( target, previousNode ) ) )
		{
			times[target] = nodeTime;

			if ( previousNode != twinSource )
			{
				previous[target] = previousNode;
			}
			else
			{
				previous[target] = twinTarget;
			}
		}

		previous[source] = target;
	}

	targetTime = times[twinTarget];
	nodeTime = times[twinSource];

	if ( nodeTime != -1 )
	{
		if ( targetTime == -1 || targetTime > nodeTime || ( targetTime == nodeTime && !isPreviousToNode ( twinTarget, twinSource ) ) )
		{
			times[twinTarget] = nodeTime;
			previous[twinTarget] = twinSource;
		}
	}

	remapNodeInwardReferencesOntoNode ( twinSource, twinTarget );
}

static void remapBackOfNodeArcsOntoNeighbour ( unsigned int source, unsigned int target )
{
	ARC * arc;
	remapNodeArcsOntoTarget ( getTwinEdge ( source ), getTwinEdge ( target ) );

	for ( arc = edge_array[source].arcs; arc != NULL; arc = arc->next )
	{
		createAnalogousArc ( target, source, arc );
	}
}

/*************************************************
Function:
    remapBackOfNodeOntoNeighbour
Description:
    Splits the edge 'source' and updates the info.
Input:
    1. source:          the source edge index
    2. sourceMarker:        the linked of the source path
    3. target:          the target edge index
    4. targetMarker:        the linked of the target path
    5. slowToFast:          1 for slow edge map onto fast edge; 0 for fast edge map onto slow edge
Output:
    None.
Return:
    None.
*************************************************/
static void remapBackOfNodeOntoNeighbour ( unsigned int source, READINTERVAL * sourceMarker, unsigned int target, READINTERVAL * targetMarker, boolean slowToFast )
{
	int offset;
	offset = remapBackOfNodeMarkersOntoNeighbour ( source, sourceMarker, target, targetMarker, slowToFast );
	remapBackOfNodeDescriptorOntoNeighbour ( source, target, slowToFast, offset );
	combineCOV ( source, edge_array[source].length, target, edge_array[target].length );
	remapBackOfNodeTimesOntoNeighbour ( source, target );
	remapBackOfNodeArcsOntoNeighbour ( source, target );
	remapNodeFibHeapReferencesOntoNode ( getTwinEdge ( source ), getTwinEdge ( target ) );
	//why not "remapNodeFibHeapReferencesOntoNode(source,target);"
	//because the downstream part of source still retains, which can serve as previousNode as before

	if ( getTwinEdge ( source ) == startingNode )
	{
		startingNode = getTwinEdge ( target );
	}
}

/*************************************************
Function:
    markerLeadsToNode
Description:
    Checks if the path exists the edge.
Input:
    1. marker:      the linked of the path
    2. node:        the edge index
Output:
    None.
Return:
    True if the edge is on the path.
*************************************************/
static boolean markerLeadsToNode ( READINTERVAL * marker, unsigned int node )
{
	READINTERVAL * currentMarker;

	for ( currentMarker = marker; currentMarker != NULL; currentMarker = currentMarker->nextInRead )
		if ( currentMarker->edgeid == node )
		{
			return true;
		}

	return false;
}

/*************************************************
Function:
    reduceNode
Description:
    Deletes the edge.
Input:
    1. node:        the edge index
Output:
    None.
Return:
    None.
*************************************************/
static void reduceNode ( unsigned int node )
{
	unsigned int bal_ed = getTwinEdge ( node );
	edge_array[node].length = 0;
	edge_array[bal_ed].length = 0;
}

/*************************************************
Function:
    reduceSlowNodes
Description:
    Removes the edges from the slow path until the edge 'finish'.
Input:
    1. slowMarker:      the linked of the slow path
    2. finish:      the edge index
Output:
    None.
Return:
    None.
*************************************************/
static void reduceSlowNodes ( READINTERVAL * slowMarker, unsigned int finish )
{
	READINTERVAL * marker;

	for ( marker = slowMarker; marker->edgeid != finish; marker = marker->nextInRead )
	{
		reduceNode ( marker->edgeid );
	}
}

static boolean markerLeadsToArc ( READINTERVAL * marker, unsigned int nodeA, unsigned int nodeB )
{
	READINTERVAL * current, *next;
	unsigned int twinA = getTwinEdge ( nodeA );
	unsigned int twinB = getTwinEdge ( nodeB );
	current = marker;

	while ( current != NULL )
	{
		next = current->nextInRead;

		if ( current->edgeid == nodeA && next->edgeid == nodeB )
		{
			return true;
		}

		if ( current->edgeid == twinB && next->edgeid == twinA )
		{
			return true;
		}

		current = next;
	}

	return false;
}

static void remapEmptyPathArcsOntoMiddlePathSimple ( READINTERVAL * emptyPath, READINTERVAL * targetPath )
{
	READINTERVAL * pathMarker, *marker;
	unsigned int start = emptyPath->prevInRead->edgeid;
	unsigned int finish = emptyPath->edgeid;
	unsigned int previousNode = start;
	unsigned int currentNode;
	ARC * originalArc = getArcBetween ( start, finish );

	if ( !originalArc )
	{
		fprintf ( stderr, "RemapEmptyPathArcsOntoMiddlePathSimple: no arc between %d and %d.\n", start, finish );
		marker = fastPath;
		fprintf ( stderr, "Fast path: " );

		while ( marker )
		{
			fprintf ( stderr, "%d,", marker->edgeid );
			marker = marker->nextInRead;
		}

		fprintf ( stderr, "\n" );
		marker = slowPath;
		fprintf ( stderr, "Slow path: " );

		while ( marker )
		{
			fprintf ( stderr, "%d,", marker->edgeid );
			marker = marker->nextInRead;
		}

		fprintf ( stderr, "\n" );
	}

	for ( pathMarker = targetPath; pathMarker->edgeid != finish; pathMarker = pathMarker->nextInRead )
	{
		currentNode = pathMarker->edgeid;
		createAnalogousArc ( previousNode, currentNode, originalArc );
		previousNode = currentNode;
	}

	createAnalogousArc ( previousNode, finish, originalArc );
	destroyArc ( start, originalArc );
}

static void remapEmptyPathMarkersOntoMiddlePathSimple ( READINTERVAL * emptyPath, READINTERVAL * targetPath, boolean slowToFast )
{
	READINTERVAL * marker, *newMarker, *previousMarker, *pathMarker, *bal_marker;
	unsigned int start = emptyPath->prevInRead->edgeid;
	unsigned int finish = emptyPath->edgeid;
	unsigned int markerStart, bal_ed;
	READINTERVAL * oldMarker = edge_array[finish].rv;

	while ( oldMarker )
	{
		marker = oldMarker;
		oldMarker = marker->nextOnEdge;
		newMarker = marker->prevInRead;

		if ( newMarker->edgeid != start )
		{
			continue;
		}

		if ( ( slowToFast && marker->readid != 2 ) || ( !slowToFast && marker->readid != 1 ) )
		{
			continue;
		}

		markerStart = marker->start;

		for ( pathMarker = targetPath; pathMarker->edgeid != finish; pathMarker = pathMarker->nextInRead )
		{
			previousMarker = newMarker;
			//maker a new marker
			newMarker = allocateRV ( marker->readid, pathMarker->edgeid );
			newMarker->start = markerStart;
			edge_array[pathMarker->edgeid].rv = addRv ( edge_array[pathMarker->edgeid].rv, newMarker );
			//maker the twin marker
			bal_ed = getTwinEdge ( pathMarker->edgeid );
			bal_marker = allocateRV ( -marker->readid, bal_ed );
			bal_marker->start = markerStart;
			edge_array[bal_ed].rv = addRv ( edge_array[bal_ed].rv, bal_marker );
			newMarker->bal_rv = bal_marker;
			bal_marker->bal_rv = newMarker;
			connectInRead ( previousMarker, newMarker );
		}

		connectInRead ( newMarker, marker );
	}
}

static void remapNodeTimesOntoForwardMiddlePath ( unsigned int source, READINTERVAL * path )
{
	READINTERVAL * marker;
	unsigned int target;
	Time nodeTime = times[source];
	unsigned int previousNode = previous[source];
	Time targetTime;

	for ( marker = path; marker->edgeid != source; marker = marker->nextInRead )
	{
		target = marker->edgeid;
		targetTime = times[target];

		if ( targetTime == -1 || targetTime > nodeTime || ( targetTime == nodeTime && !isPreviousToNode ( target, previousNode ) ) )
		{
			times[target] = nodeTime;
			previous[target] = previousNode;
		}

		previousNode = target;
	}

	previous[source] = previousNode;
}

static void remapNodeTimesOntoTwinMiddlePath ( unsigned int source, READINTERVAL * path )
{
	READINTERVAL * marker;
	unsigned int target;
	unsigned int previousNode = getTwinEdge ( source );
	Time targetTime;
	READINTERVAL * limit = path->prevInRead->bal_rv;
	Time nodeTime = times[limit->edgeid];
	marker = path;

	while ( marker->edgeid != source )
	{
		marker = marker->nextInRead;
	}

	marker = marker->bal_rv;

	while ( marker != limit )
	{
		marker = marker->nextInRead;
		target = marker->edgeid;
		targetTime = times[target];

		if ( targetTime == -1 || targetTime > nodeTime || ( targetTime == nodeTime && !isPreviousToNode ( target, previousNode ) ) )
		{
			times[target] = nodeTime;
			previous[target] = previousNode;
		}

		previousNode = target;
	}
}

/*************************************************
Function:
    remapEmptyPathOntoMiddlePath
Description:
    Updates the info of the slow path according to the fast path.
Input:
    1. emptyPath:       the linked of the path
    2. targetPath:      the linked of the path
    3. slowToFast:      0 for 'emptyPath' is the fast path and 'targetPath' is the slow path
Output:
    None.
Return:
    None.
*************************************************/
static void remapEmptyPathOntoMiddlePath ( READINTERVAL * emptyPath, READINTERVAL * targetPath, boolean slowToFast )
{
	unsigned int start = emptyPath->prevInRead->edgeid;
	unsigned int finish = emptyPath->edgeid;

	// Remapping markers
	if ( !markerLeadsToArc ( targetPath, start, finish ) )
	{
		remapEmptyPathArcsOntoMiddlePathSimple ( emptyPath, targetPath );
	}

	remapEmptyPathMarkersOntoMiddlePathSimple ( emptyPath, targetPath, slowToFast );

	//Remap times and previous(if necessary)
	if ( getNodePrevious ( finish ) == start )
	{
		remapNodeTimesOntoForwardMiddlePath ( finish, targetPath );
	}

	if ( getNodePrevious ( getTwinEdge ( start ) ) == getTwinEdge ( finish ) )
	{
		remapNodeTimesOntoTwinMiddlePath ( finish, targetPath );
	}
}

/*************************************************
Function:
    cleanUpRedundancy
Description:
    Merges the bubble and copies the position info from the slow path to the fast path.
Input:
    None.
Output:
    None.
Return:
    1 if the bubble is merged successfully.
*************************************************/
static boolean cleanUpRedundancy ()
{
	READINTERVAL * slowMarker = slowPath->nextInRead, *fastMarker = fastPath->nextInRead;
	unsigned int slowNode, fastNode;
	int slowLength, fastLength;
	int fastConstraint = 0;
	int slowConstraint = 0;
	int finalLength;
	attachPath ( slowPath );
	attachPath ( fastPath );
	mapSlowOntoFast ();
	finalLength = mapDistancesOntoPaths ();
	slowLength = fastLength = 0;

	while ( slowMarker != NULL && fastMarker != NULL )
	{
		if ( !slowMarker->nextInRead )
		{
			slowLength = finalLength;
		}
		else
		{
			slowLength = slowToFastMapping[slowMarker->bal_rv->start - 1];

			if ( slowLength < slowConstraint )
			{
				slowLength = slowConstraint;
			}
		}

		fastLength = fastMarker->bal_rv->start - 1;

		if ( fastLength < fastConstraint )
		{
			fastLength = fastConstraint;
		}

		slowNode = slowMarker->edgeid;
		fastNode = fastMarker->edgeid;

		if ( false )
			{ fprintf ( stderr, "Slow %d    Fast %d.\n", slowLength, fastLength ); }

		if ( slowNode == fastNode )
		{
			if ( false )
				{ fprintf ( stderr, "0/ Already merged together %d == %d.\n", slowNode, fastNode ); }

			if ( fastLength > slowLength )
			{
				slowConstraint = fastLength;
			}

			fastConstraint = slowLength;
			slowMarker = slowMarker->nextInRead;
			fastMarker = fastMarker->nextInRead;
		}
		else if ( slowNode == getTwinEdge ( fastNode ) )
		{
			if ( false )
				{ fprintf ( stderr, "1/ Creme de la hairpin %d $$ %d.\n", slowNode, fastNode ); }

			if ( fastLength > slowLength )
			{
				slowConstraint = fastLength;
			}

			fastConstraint = slowLength;
			slowMarker = slowMarker->nextInRead;
			fastMarker = fastMarker->nextInRead;
		}
		else if ( markerLeadsToNode ( slowMarker, fastNode ) )
		{
			if ( false )
			{
				fprintf ( stderr, "2/ Remapping empty fast arc onto slow nodes.\n" );
			}

			reduceSlowNodes ( slowMarker, fastNode );
			remapEmptyPathOntoMiddlePath ( fastMarker, slowMarker, FAST_TO_SLOW );

			while ( slowMarker->edgeid != fastNode )
			{
				slowMarker = slowMarker->nextInRead;
			}
		}
		else if ( markerLeadsToNode ( fastMarker, slowNode ) )
		{
			if ( false )
			{
				fprintf ( stderr, "3/ Remapping empty slow arc onto fast nodes.\n" );
			}

			remapEmptyPathOntoMiddlePath ( slowMarker, fastMarker, SLOW_TO_FAST );

			while ( fastMarker->edgeid != slowNode )
			{
				fastMarker = fastMarker->nextInRead;
			}
		}
		else if ( slowLength == fastLength )
		{
			if ( false )
			{
				fprintf ( stderr, "A/ Mapped equivalent nodes together %d <=> %d.\n", slowNode, fastNode );
			}

			remapNodeOntoNeighbour ( slowNode, fastNode );
			slowMarker = slowMarker->nextInRead;
			fastMarker = fastMarker->nextInRead;
		}
		else if ( slowLength < fastLength )
		{
			if ( false )
			{
				fprintf ( stderr, "B/ Mapped back of fast node into slow %d -> %d.\n", fastNode, slowNode );
			}

			remapBackOfNodeOntoNeighbour ( fastNode, fastMarker, slowNode, slowMarker, FAST_TO_SLOW );
			slowMarker = slowMarker->nextInRead;
		}
		else
		{
			if ( false )
			{
				fprintf ( stderr, "C/ Mapped back of slow node into fast %d -> %d.\n", slowNode, fastNode );
			}

			remapBackOfNodeOntoNeighbour ( slowNode, slowMarker, fastNode, fastMarker, SLOW_TO_FAST );
			fastMarker = fastMarker->nextInRead;
		}
	}

	detachPath ( fastPath );
	detachPath ( slowPath );
	return 1;
}

/*************************************************
Function:
    comparePaths
Description:
    1. Gets two path of the bubble.
    2. Checks if the bubble is suitable to merge.
    3. Merges the bubble.
Input:
    1. destination:     the edge index, the start node of the fast path
    2. origin:      the edge index, the start node of the slow path
Output:
    None.
Return:
    None.
*************************************************/
static void comparePaths ( unsigned int destination, unsigned int origin )
{
	int slowLength, fastLength, i;
	unsigned int fastNode, slowNode;
	READINTERVAL * marker;
	slowLength = fastLength = 0;
	fastNode = destination;
	slowNode = origin;
	btCounter++;

	while ( fastNode != slowNode )
	{
		if ( times[fastNode] > times[slowNode] )
		{
			fastLength++;
			fastNode = previous[fastNode];
		}
		else if ( times[fastNode] < times[slowNode] )
		{
			slowLength++;
			slowNode = previous[slowNode];
		}
		else if ( isPreviousToNode ( slowNode, fastNode ) )
		{
			while ( fastNode != slowNode )
			{
				fastLength++;
				fastNode = previous[fastNode];
			}
		}
		else if ( isPreviousToNode ( fastNode, slowNode ) )
		{
			while ( slowNode != fastNode )
			{
				slowLength++;
				slowNode = previous[slowNode];
			}
		}
		else
		{
			fastLength++;
			fastNode = previous[fastNode];
			slowLength++;
			slowNode = previous[slowNode];
		}

		if ( slowLength > MAXNODELENGTH || fastLength > MAXNODELENGTH )
		{
			return;
		}
	}

	if ( fastLength == 0 )
	{
		return;
	}

	marker = allocateRV ( 1, destination );
	fastPath = marker;

	for ( i = 0; i < fastLength; i++ )
	{
		marker = allocateRV ( 1, previous[fastPath->edgeid] );
		marker->nextInRead = fastPath;
		fastPath->prevInRead = marker;
		fastPath = marker;
	}

	marker = allocateRV ( 2, destination );
	slowPath = marker;
	marker = allocateRV ( 2, origin );
	marker->nextInRead = slowPath;
	slowPath->prevInRead = marker;
	slowPath = marker;

	for ( i = 0; i < slowLength; i++ )
	{
		marker = allocateRV ( 2, previous[slowPath->edgeid] );
		marker->nextInRead = slowPath;
		slowPath->prevInRead = marker;
		slowPath = marker;
	}

	fastSeqLength = extractSequence ( fastPath, fastSequence );
	slowSeqLength = extractSequence ( slowPath, slowSequence );

	/*
	   if(destination==6359){
	   printf("destination %d, slowLength %d, fastLength %d\n",destination,slowLength,fastLength);
	   printf("fastSeqLength %d, slowSeqLength %d\n",fastSeqLength,slowSeqLength);
	   }
	 */
	if ( !fastSeqLength || !slowSeqLength )
	{
		detachPathSingle ( slowPath );
		detachPathSingle ( fastPath );
		return;
	}

	cmpCounter++;

	if ( !compareSequences ( fastSequence, slowSequence, fastSeqLength, slowSeqLength ) )
	{
		detachPathSingle ( slowPath );
		detachPathSingle ( fastPath );
		return;
	}

	//only merge clean bubble ...
	if ( clean )
	{
		unsigned int bal_ed;
		unsigned int arcRight_n, arcLeft_n;
		READINTERVAL * tmp;
		tmp = fastPath->nextInRead;

		while ( tmp->nextInRead )
		{
			bal_ed = getTwinEdge ( tmp->edgeid );
			arcCount ( tmp->edgeid, &arcRight_n );
			arcCount ( bal_ed, &arcLeft_n );

			if ( arcRight_n != 1 || arcLeft_n != 1 ) //not clean bubble
			{
				return;
			}

			tmp = tmp->nextInRead;
		}

		tmp = slowPath->nextInRead;

		while ( tmp->nextInRead )
		{
			bal_ed = getTwinEdge ( tmp->edgeid );
			arcCount ( tmp->edgeid, &arcRight_n );
			arcCount ( bal_ed, &arcLeft_n );

			if ( arcRight_n != 1 || arcLeft_n != 1 ) //not clean bubble
			{
				return;
			}

			tmp = tmp->nextInRead;
		}
	}

	simiCounter++;
	pinCounter += cleanUpRedundancy ();

	if ( pinCounter % 100000 == 0 )
	{
		fprintf ( stderr, ".............%lld bubbles merged.\n", pinCounter );
	}

	HasChanged = 1;
}

/*************************************************
Function:
    tourBusArc
Description:
    1. Searchs the bubble using  the "Tour Bus algorithm" by the Function "tourBus", "tourBusNode", "tourBusArc".
    2. Checks if the edge "destination" is visited. if not, return.
    3. Compares two path if the edge "destination" is visited twice.
Input:
    1. origin:          the origin edge
    2. destination: the destination edge
    3. arcMulti:        the weight of the arc, which is from the origin edge to the destination edge
    4. originTime:      the times to the origin edge
Output:
    None.
Return:
    None.
*************************************************/

static void tourBusArc ( unsigned int origin, unsigned int destination, unsigned int arcMulti, Time originTime )
{
	Time arcTime, totalTime, destinationTime;
	unsigned int oldPrevious = previous[destination];

	if ( oldPrevious == origin || edge_array[destination].multi == 1 )
	{
		return;
	}

	arcCounter++;

	if ( arcMulti > 0 )
	{
		arcTime = ( ( Time ) edge_array[origin].length ) / ( ( Time ) arcMulti );
	}
	else
	{
		arcTime = 0.0;
		fprintf ( stderr, "Arc from %d to %d with flags %d originTime %f, arc %d.\n", origin, destination, edge_array[destination].multi, originTime, arcMulti );
	}

	totalTime = originTime + arcTime;
	/*
	   if(destination==289129||destination==359610){
	   printf("arc from %d to %d with flags %d time %f originTime %f, arc %d\n",
	   origin,destination,edge_array[destination].multi,totalTime,originTime,arcMulti);
	   fflush(stdout);
	   }
	 */
	destinationTime = times[destination];

	if ( destinationTime == -1 )
	{
		times[destination] = totalTime;
		dheapNodes[destination] = insertNodeIntoDHeap ( dheap, totalTime, destination );
		dnodeCounter++;
		previous[destination] = origin;
		return;
	}
	else if ( destinationTime > totalTime )
	{
		if ( dheapNodes[destination] == NULL )
		{
			return;
		}

		replaceCounter++;
		times[destination] = totalTime;
		replaceKeyInDHeap ( dheap, dheapNodes[destination], totalTime );
		previous[destination] = origin;
		comparePaths ( destination, oldPrevious );
		return;
	}
	else
	{
		if ( destinationTime == times[origin] && isPreviousToNode ( destination, origin ) )
		{
			return;
		}

		comparePaths ( destination, origin );
	}
}

/*************************************************
Function:
    tourBusNode
Description:
    Searchs the bubble using the "Tour Bus algorithm" by the function "tourBus", "tourBusNode", "tourBusArc".
Input:
    1. node:        the current edge
Output:
    None.
Return:
    None.
*************************************************/
static void tourBusNode ( unsigned int node )
{
	ARC * parc;
	int index = 0, outNodeNum;
	expanded[expCounter++] = node;
	activeNode = node;
	parc = edge_array[activeNode].arcs;

	while ( parc )
	{
		outArcArray[index] = parc;
		outNodeArray[index++] = parc->to_ed;

		if ( index >= MAXCONNECTION )
		{
			break;
		}

		parc = parc->next;
	}

	outNodeNum = index;
	HasChanged = 0;

	for ( index = 0; index < outNodeNum; index++ )
	{
		if ( HasChanged )
		{
			parc = getArcBetween ( activeNode, outNodeArray[index] );
			getArcCounter++;
		}
		else
		{
			parc = outArcArray[index];
		}

		if ( !parc )
		{
			continue;
		}

		tourBusArc ( activeNode, outNodeArray[index], parc->multiplicity, times[activeNode] );
	}
}

/*
static void dumpNodeFromDHeap()
{
    unsigned int currentNode;

    while((currentNode = removeNextNodeFromDHeap(dheap))!=0){
        rnodeCounter++;
        times[currentNode] = -1;
        previous[currentNode] = 0;
        dheapNodes[currentNode] = NULL;
        if(dnodeCounter-rnodeCounter<250)
            break;
    }
}
*/

/*************************************************
Function:
    tourBus
Description:
    Searchs the bubble using the "Tour Bus algorithm" by the function "tourBus", "tourBusNode", "tourBusArc".
Input:
    1. startingPoint:       the start edge
Output:
    None.
Return:
    None.
*************************************************/
static void tourBus ( unsigned int startingPoint )
{
	unsigned int currentNode = startingPoint;
	times[startingPoint] = 0;
	previous[startingPoint] = currentNode;

	while ( currentNode > 0 )
	{
		dheapNodes[currentNode] = NULL;
		tourBusNode ( currentNode );
		currentNode = removeNextNodeFromDHeap ( dheap );

		if ( currentNode > 0 )
		{
			rnodeCounter++;
		}
	}
}

/*************************************************
Function:
    bubblePinch
Description:
    Removes bubbles with the Tour Bus algorithm.
Input:
    1. simiCutoff:      the minimum requirements for similarity
    2. outfile:     the output prefix
    3. M:           the strength of merging bubble
    4. isIter:      whether it's multikmer
    5. last:        whether it's the last iteration
Output:
    None.
Return:
    None.
*************************************************/
void bubblePinch ( double simiCutoff, char * outfile, int M, boolean isIter, boolean last )
{
	new_num_vt = 2 * num_vt;
	unsigned int index, counter = 0;
	unsigned int startingNode;
	char temp[256];
	sprintf ( temp, "%s.pathpair", outfile );
	caseA = caseB = caseC = caseD = caseE = 0;
	progress = 0;
	arcCounter = 0;
	dnodeCounter = 0;
	rnodeCounter = 0;
	btCounter = 0;
	cmpCounter = 0;
	simiCounter = 0;
	pinCounter = 0;
	replaceCounter = 0;
	getArcCounter = 0;
	cutoff = 1.0 - simiCutoff;

	if ( M <= 1 )
	{
		MAXNODELENGTH = 3;
		DIFF = 2;
	}
	else if ( M == 2 )
	{
		MAXNODELENGTH = 9;
		DIFF = 3;
	}
	else
	{
		MAXNODELENGTH = 30;
		DIFF = 10;
	}

	fprintf ( stderr, "Start to pinch bubbles, cutoff %f, MAX NODE NUM %d, MAX DIFF NUM %d.\n", cutoff, MAXNODELENGTH, DIFF );
	createRVmemo ();
	times = ( Time * ) ckalloc ( ( num_ed + 1 ) * sizeof ( Time ) );
	previous = ( unsigned int * ) ckalloc ( ( num_ed + 1 ) * sizeof ( unsigned int ) );
	expanded = ( unsigned int * ) ckalloc ( ( num_ed + 1 ) * sizeof ( unsigned int ) );
	dheapNodes = ( DFibHeapNode ** ) ckalloc ( ( num_ed + 1 ) * sizeof ( DFibHeapNode * ) );
	WORDFILTER = createFilter ( overlaplen );

	for ( index = 1; index <= num_ed; index++ )
	{
		times[index] = -1;
		previous[index] = 0;
		dheapNodes[index] = NULL;
	}

	dheap = newDFibHeap ();
	eligibleStartingPoints = ( unsigned int * ) ckalloc ( ( num_ed + 1 ) * sizeof ( unsigned int ) );
	resetNodeStatus ();
	createArcLookupTable ();
	recordArcsInLookupTable ();

	while ( ( startingNode = nextStartingPoint () ) > 0 )
	{
		counter++;
		expCounter = 0;
		tourBus ( startingNode );
		updateNodeStatus ();
	}

	resetNodeStatus ();
	deleteArcLookupTable ();
	destroyReadIntervMem ();
	fprintf ( stderr, "%d start points, %lld dheap nodes.\n", counter, dnodeCounter );
	fprintf ( stderr, "%lld pair(s) found, %lld pair of path(s) compared, %lld pair(s) merged.\n", btCounter, cmpCounter, pinCounter );
	fprintf ( stderr, "Sequence comparison failed:\n" );
	fprintf ( stderr, " Path crossing deleted edge                         %lld\n", caseA );
	fprintf ( stderr, " Length difference of two paths greater than two    %lld\n", caseB );
	fprintf ( stderr, " Mismatch score greater than cutoff (%d)             %lld\n", DIFF, caseC );
	fprintf ( stderr, " Mismatch score ratio greater than cutoff (%.1f)     %lld\n", cutoff, caseD );
	fprintf ( stderr, " Path length shorter than (Kmer-1)                  %lld\n", caseE );
	free ( ( void * ) eligibleStartingPoints );
	destroyDHeap ( dheap );
	free ( ( void * ) dheapNodes );
	free ( ( void * ) times );
	free ( ( void * ) previous );
	free ( ( void * ) expanded );
	linearConcatenate ( isIter, last );
}

/*
 * concatenateEdge.c
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

void copySeq ( char * targetS, char * sourceS, int pos, int length )
{
	char ch;
	int i, index;
	index = pos;

	for ( i = 0; i < length; i++ )
	{
		ch = getCharInTightString ( sourceS, i );
		writeChar2tightString ( ch, targetS, index++ );
	}
}

/*************************************************
Function:
    linearUpdateConnection
Description:
    A path from e1 to e2 is merged in to e1 (indicate=0) or e2 (indicate=1), update graph topology.
Input:
    1. e1:      edge1 index
    2. e2:      edge2 index
    3. indicate:    indicate which edge to be merged into
Output:
    None.
Return:
    None.
*************************************************/
void linearUpdateConnection ( unsigned int e1, unsigned int e2, int indicate )
{
	unsigned int bal_ed;
	ARC * parc;

	if ( !indicate )
	{
		edge_array[e1].to_vt = edge_array[e2].to_vt;
		bal_ed = getTwinEdge ( e1 );
		parc = edge_array[e2].arcs;

		while ( parc )
		{
			parc->bal_arc->to_ed = bal_ed;
			parc = parc->next;
		}

		edge_array[e1].arcs = edge_array[e2].arcs;
		edge_array[e2].arcs = NULL;

		if ( edge_array[e1].length || edge_array[e2].length )
			{ edge_array[e1].cvg = ( edge_array[e1].cvg * edge_array[e1].length + edge_array[e2].cvg * edge_array[e2].length ) / ( edge_array[e1].length + edge_array[e2].length ); }

		edge_array[e2].deleted = 1;
	}
	else
	{
		//all the arcs pointing to e1 switch to e2
		parc = edge_array[getTwinEdge ( e1 )].arcs;

		while ( parc )
		{
			parc->bal_arc->to_ed = e2;
			parc = parc->next;
		}

		edge_array[e1].arcs = NULL;
		edge_array[e2].from_vt = edge_array[e1].from_vt;

		if ( edge_array[e1].length || edge_array[e2].length )
			{ edge_array[e2].cvg = ( edge_array[e1].cvg * edge_array[e1].length + edge_array[e2].cvg * edge_array[e2].length ) / ( edge_array[e1].length + edge_array[e2].length ); }

		edge_array[e1].deleted = 1;
	}
}

static void printEdgeSeq ( FILE * fp, char * tightSeq, int len )
{
	int i;

	for ( i = 0; i < len; i++ )
	{
		fprintf ( fp, "%c", int2base ( ( int ) getCharInTightString ( tightSeq, i ) ) );

		if ( ( i + overlaplen + 1 ) % 100 == 0 && i < len - 1 )
		{
			fprintf ( fp, "\n" );
		}
	}

	fprintf ( fp, "\n" );
}
void allpathUpdateEdge ( unsigned int e1, unsigned int e2, int indicate, boolean last )
{
	int tightLen;
	char * tightSeq = NULL;

	if ( edge_array[e1].cvg == 0 )
	{
		edge_array[e1].cvg = edge_array[e2].cvg;
	}

	if ( edge_array[e2].cvg == 0 )
	{
		edge_array[e2].cvg = edge_array[e1].cvg;
	}

	/*
	   if(edge_array[e1].length&&edge_array[e2].length){
	   fprintf(stderr,">e1\n");
	   printEdgeSeq(stderr,edge_array[e1].seq,edge_array[e1].length);
	   fprintf(stderr,">e2\n");
	   printEdgeSeq(stderr,edge_array[e2].seq,edge_array[e2].length);
	   } */
	unsigned int cvgsum = edge_array[e1].cvg * edge_array[e1].length + edge_array[e2].cvg * edge_array[e2].length;
	tightLen = edge_array[e1].length + edge_array[e2].length;

	if ( tightLen )
	{
		tightSeq = ( char * ) ckalloc ( ( tightLen / 4 + 1 ) * sizeof ( char ) );
	}

	tightLen = 0;

	if ( edge_array[e1].length )
	{
		copySeq ( tightSeq, edge_array[e1].seq, 0, edge_array[e1].length );
		tightLen = edge_array[e1].length;

		if ( edge_array[e1].seq )
		{
			free ( ( void * ) edge_array[e1].seq );
			edge_array[e1].seq = NULL;
		}
		else
		{
			fprintf ( stderr, "AllpathUpdateEdge: edge %d with length %d, but without seq.\n", e1, edge_array[e1].length );
		}
	}

	if ( edge_array[e2].length )
	{
		copySeq ( tightSeq, edge_array[e2].seq, tightLen, edge_array[e2].length );
		tightLen += edge_array[e2].length;

		if ( edge_array[e2].seq )
		{
			free ( ( void * ) edge_array[e2].seq );
			edge_array[e2].seq = NULL;
		}
		else
		{
			fprintf ( stderr, "AllpathUpdateEdge: edge %d with length %d, but without seq.\n", e2, edge_array[e2].length );
		}
	}

	/*
	   if(edge_array[e1].length&&edge_array[e2].length){
	   fprintf(stderr,">e1+e2\n");
	   printEdgeSeq(stderr,tightSeq,tightLen);
	   }
	 */
	//edge_array[e2].extend_len = tightLen-edge_array[e2].length;
	//the sequence of e1 is to be updated
	if ( !indicate )
	{
		edge_array[e2].length = 0;  //e1 is removed from the graph
		edge_array[e1].to_vt = edge_array[e2].to_vt;    //e2 is part of e1 now
		edge_array[e1].length = tightLen;
		edge_array[e1].seq = tightSeq;

		if ( tightLen )
		{
			edge_array[e1].cvg = cvgsum / tightLen;
		}

		if ( last )
			{ edge_array[e1].cvg = edge_array[e1].cvg > 0 ? edge_array[e1].cvg : 1; }

		//      edge_array[e1].cvg = edge_array[e1].cvg > 0 ? edge_array[e1].cvg : 1;
	}
	else
	{
		edge_array[e1].length = 0;  //e1 is removed from the graph
		edge_array[e2].from_vt = edge_array[e1].from_vt;    //e1 is part of e2 now
		edge_array[e2].length = tightLen;
		edge_array[e2].seq = tightSeq;

		if ( tightLen )
		{
			edge_array[e2].cvg = cvgsum / tightLen;
		}

		if ( last )
			{ edge_array[e2].cvg = edge_array[e2].cvg > 0 ? edge_array[e2].cvg : 1; }

		//      edge_array[e2].cvg = edge_array[e2].cvg > 0 ? edge_array[e2].cvg : 1;
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
    linearConcatenate
Description:
    Concatenates two edges if they are linearly linked.
Input:
    1. isIter:      whether it is multi kmer
    2. last:        whether it is the last iteration
Output:
    None.
Return:
    None.
*************************************************/
void linearConcatenate ( boolean isIter, boolean last )
{
	unsigned int i;
	int conc_c = 1;
	int counter;
	unsigned int from_ed, to_ed, bal_ed;
	ARC * parc, *parc2;
	unsigned int bal_fe;
	int donot1 = 0;
	int round = 1;

	//debugging(30514);
	while ( conc_c )
	{
		conc_c = 0;
		counter = 0;
		donot1 = 0;

		for ( i = 1; i <= num_ed; i++ ) //num_ed
		{
			if ( edge_array[i].deleted || EdSameAsTwin ( i ) )
			{
				continue;
			}

			if ( edge_array[i].length > 0 )
			{
				counter++;
			}

			parc = edge_array[i].arcs;

			if ( !parc || parc->next )
			{
				continue;
			}

			to_ed = parc->to_ed;
			bal_ed = getTwinEdge ( to_ed );
			parc2 = edge_array[bal_ed].arcs;

			if ( bal_ed == to_ed || !parc2 || parc2->next )
			{
				continue;
			}

			from_ed = i;

			if ( from_ed == to_ed || from_ed == bal_ed )
			{
				continue;
			}

			if ( parc->multiplicity <= arcfilter )
			{
				donot1++;
				continue;
			}

			//linear connection found
			conc_c++;
			linearUpdateConnection ( from_ed, to_ed, 0 );
			allpathUpdateEdge ( from_ed, to_ed, 0, last );
			bal_fe = getTwinEdge ( from_ed );
			linearUpdateConnection ( bal_ed, bal_fe, 1 );
			allpathUpdateEdge ( bal_ed, bal_fe, 1, last );
			/*
			   if(from_ed==6589||to_ed==6589)
			   printf("%d <- %d (%d)\n",from_ed,to_ed,i);
			   if(bal_fe==6589||bal_ed==6589)
			   printf("%d <- %d (%d)\n",bal_fe,bal_ed,i);
			 */
		}

		fprintf ( stderr, "%d edge(s) concatenated in cycle %d.\n", conc_c, round++ );

		if ( arcfilter )
			{ fprintf ( stderr, "%d edge(s) with weight %d were not linearized.\n", donot1, arcfilter ); }
	}
}

/*
 * linearEdge.c
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

#ifdef MER127
//Get the char in Kmer.
char getCharInKmer ( Kmer kmer, int pos )
{
	if ( 2 * pos < 64 )
	{
		kmer.low2 >>= 2 * pos;
		return kmer.low2 & 3;
	}
	else if ( 2 * pos < 128 )
	{
		kmer.high2 >>= 2 * pos - 64;
		return kmer.high2 & 3;
	}
	else if ( 2 * pos < 192 )
	{
		kmer.low1 >>= 2 * pos - 128;
		return kmer.low1 & 3;
	}
	else
	{
		kmer.high1 >>= 2 * pos - 192;
		return kmer.high1 & 3;
	}
}
#else
char getCharInKmer ( Kmer kmer, int pos )
{
	if ( 2 * pos < 64 )
	{
		kmer.low >>= 2 * pos;
		return kmer.low & 3;
	}
	else
	{
		kmer.high >>= 2 * pos - 64;
		return kmer.high & 3;
	}
}
#endif

/*************************************************
Function:
    copyinter
Description:
    Coyies the interger kmer to char seq.
Input:
    1. targetS:     the target seq
    2. sourceS:     the source kmer
    3. pos:         the pos at kmer
    4. length:      the length of sequence to copy
Output:
    None.
Return:
    None.
*************************************************/
void copyinter ( char * targetS, Kmer sourceS, int pos, int length )
{
	char ch;
	int i, index;
	index = pos;

	for ( i = 0; i < length; ++i )
	{
		ch = getCharInKmer ( sourceS, step - i - 1 );
		writeChar2tightString ( ch, targetS, index++ );
	}
}

/*************************************************
Function:
    copySeq2
Description:
    Copies seq from source to target.
Input:
    1. targetS:     the target seq
    2. sourceS:     the source seq
    3. pos:         the pos at source
    4. length:      the length of sequence to copy
Output:
    None.
Return:
    None.
*************************************************/
void copySeq2 ( char * targetS, char * sourceS, int pos, int length )
{
	char ch;
	int i, index;
	index = pos;

	for ( i = 0; i < length; ++i )
	{
		ch = getCharInTightString ( sourceS, i );
		writeChar2tightString ( ch, targetS, index++ );
	}
}

/*************************************************
Function:
    checkstep
Description:
    Checks the step between two kmer.
Input:
    1. to_vt:           the to vertex
    2. from_vt:     the from vertex
Output:
    None.
Return:
    step between two kmer.
*************************************************/
int checkstep ( unsigned int to_vt, unsigned int from_vt )
{
	Kmer to, from;
	Kmer filtertemp;
	to = vt_arraynew[to_vt].kmer;
	from = vt_arraynew[from_vt].kmer;
	int i = 1;
	filtertemp = createFilter ( overlaplen - i );

	if ( KmerEqual ( KmerRightBitMove ( from, i << 1 ), KmerAnd ( to, filtertemp ) ) )
		{ return i; }

	fprintf ( stderr, "When checking step of two edge, step is not found and step is changed to %d, 'to kmer' is ", step );
	printKmerSeq ( stderr, to );
	fprintf ( stderr, " , 'from kmer' is " );
	printKmerSeq ( stderr, from );
	fprintf ( stderr, " .\n" );
	return step;
}

/*************************************************
Function:
    linearUpdateConnection2
Description:
    Updates the arc and coverage information of edges to be merged.
Input:
    1. e1:          edge1 index
    2. e2:          edge2 index
    3. indicate:        indicates which edge would remain
Output:
    None.
Return:
    None.
*************************************************/
void linearUpdateConnection2 ( unsigned int e1, unsigned int e2, int indicate )
{
	unsigned int bal_ed;
	ARC * parc;

	//caution: length and seq
	if ( !indicate )
	{
		//      edge_array[e1].to_vt = edge_array[e2].to_vt;
		bal_ed = getTwinEdge ( e1 );
		parc = edge_array[e2].arcs;

		while ( parc )
		{
			parc->bal_arc->to_ed = bal_ed;
			parc = parc->next;
		}

		edge_array[e1].arcs = edge_array[e2].arcs;
		edge_array[e2].arcs = NULL;

		if ( edge_array[e1].length < 0 || edge_array[e2].length < 0 )
			{ fprintf ( stderr, "Error: length < 0.\n" ); }

		if ( ( edge_array[e1].length || edge_array[e2].length ) && ( edge_array[e1].length + edge_array[e2].length + step > 0 ) )
			edge_array[e1].cvg = ( edge_array[e1].cvg * ( edge_array[e1].length + step )
			                       + edge_array[e2].cvg * ( edge_array[e2].length + step ) )
			                     / ( edge_array[e1].length + edge_array[e2].length + step );

		edge_array[e2].deleted = 1;
	}
	else
	{
		//all the arcs pointing to e1 switched to e2
		parc = edge_array[getTwinEdge ( e1 )].arcs;

		while ( parc )
		{
			parc->bal_arc->to_ed = e2;
			parc = parc->next;
		}

		edge_array[e1].arcs = NULL;

		//      edge_array[e2].from_vt = edge_array[e1].from_vt;
		if ( edge_array[e1].length < 0 || edge_array[e2].length < 0 )
			{ fprintf ( stderr, "Error: length < 0.\n" ); }

		if ( ( edge_array[e1].length || edge_array[e2].length ) && ( edge_array[e1].length + edge_array[e2].length + step > 0 ) )
			edge_array[e2].cvg = ( edge_array[e1].cvg * ( edge_array[e1].length + step )
			                       + edge_array[e2].cvg * ( edge_array[e2].length + step ) )
			                     / ( edge_array[e1].length + edge_array[e2].length + step );

		edge_array[e1].deleted = 1;
	}
}

/*************************************************
Function:
    allpathUpdateEdge2
Description:
    Update graph topology.
Input:
    1. e1:          edge1 index
    2. e2:          edge2 index
    3. indicate:        indicates which edge would remain.
    4. last:            whether it's the last iterate.
Output:
    None.
Return:
    None.
*************************************************/
void allpathUpdateEdge2 ( unsigned int e1, unsigned int e2, int indicate, boolean last )
{
	int tightLen;
	char * tightSeq = NULL;
	int tempstep =  0;

	//caution: length and seq
	if ( edge_array[e1].cvg == 0 )
		{ edge_array[e1].cvg = edge_array[e2].cvg; }

	if ( edge_array[e2].cvg == 0 )
		{ edge_array[e2].cvg = edge_array[e1].cvg; }

	unsigned int cvgsum =
	    edge_array[e1].cvg * ( edge_array[e1].length + step )
	    + edge_array[e2].cvg * ( edge_array[e2].length + step );
	tightLen = edge_array[e1].length + edge_array[e2].length + step;

	if ( tightLen )
		{ tightSeq = ( char * ) ckalloc ( ( tightLen / 4 + 1 ) * sizeof ( char ) ); }

	tightLen = 0;

	if ( edge_array[e1].length )
	{
		copySeq2 ( tightSeq, edge_array[e1].seq, 0, edge_array[e1].length );
		tightLen = edge_array[e1].length;

		if ( edge_array[e1].seq )
		{
			free ( ( void * ) edge_array[e1].seq );
			edge_array[e1].seq = NULL;
		}
		else
			{ fprintf ( stderr, "AllpathUpdateEdge: edge %d with length %d, but without seq.\n", e1, edge_array[e1].length ); }
	}

	{
		if ( step > 0 )
		{
			tempstep = checkstep ( edge_array[e1].to_vt, edge_array[e2].from_vt );
			copyinter ( tightSeq, vt_arraynew[edge_array[e2].from_vt].kmer, tightLen, tempstep );
			tightLen += tempstep;
		}
	}

	if ( edge_array[e2].length )
	{
		copySeq2 ( tightSeq, edge_array[e2].seq, tightLen, edge_array[e2].length );
		tightLen += edge_array[e2].length;

		if ( edge_array[e2].seq )
		{
			free ( ( void * ) edge_array[e2].seq );
			edge_array[e2].seq = NULL;
		}
		else
			{ fprintf ( stderr, "AllpathUpdateEdge: edge %d with length %d, but without seq.\n", e2, edge_array[e2].length ); }
	}

	//edge_array[e2].extend_len = tightLen-edge_array[e2].length;
	//the sequence of e1 is to be updated
	if ( !indicate )
	{
		edge_array[e2].length = 0;    //e2 is removed from the graph
		edge_array[e1].to_vt = edge_array[e2].to_vt;      //e2 is part of e1 now
		edge_array[e1].length = tightLen;
		edge_array[e1].seq = tightSeq;

		if ( tightLen )
			{ edge_array[e1].cvg = cvgsum / tightLen; }

		if ( last )
			{ edge_array[e1].cvg = edge_array[e1].cvg > 0 ? edge_array[e1].cvg : 1; }
	}
	else
	{
		edge_array[e1].length = 0;   //e1 is removed from the graph
		edge_array[e2].from_vt = edge_array[e1].from_vt;      //e1 is part of e2 now
		edge_array[e2].length = tightLen;
		edge_array[e2].seq = tightSeq;

		if ( tightLen )
			{ edge_array[e2].cvg = cvgsum / tightLen; }

		if ( last )
			{ edge_array[e2].cvg = edge_array[e2].cvg > 0 ? edge_array[e2].cvg : 1; }
	}
}

static void debugging ( unsigned int i )
{
	ARC * parc;
	parc = edge_array[i].arcs;

	if ( !parc )
		{ fprintf ( stderr, "No downward connection for %d.\n", i ); }

	while ( parc )
	{
		fprintf ( stderr, "%d -> %d\n", i, parc->to_ed );
		parc = parc->next;
	}
}


/*************************************************
Function:
    linearConcatenate2
Description:
    Concatenates two edges if they are linearly linked.
Input:
    1. last:            whether it's the last iterations.
Output:
    None.
Return:
    None.
*************************************************/
void linearConcatenate2 ( boolean last )
{
	unsigned int i;
	int conc_c = 1;
	int counter;
	unsigned int from_ed, to_ed, bal_ed;
	ARC * parc, *parc2;
	unsigned int bal_fe;
	ARC * temp;
	int donot1 = 0;
	int round = 1;

	while ( conc_c )
	{
		conc_c = 0;
		counter = 0;
		donot1 = 0;

		for ( i = 1; i <= num_ed; i++ )
		{
			if ( edge_array[i].deleted || EdSameAsTwin ( i ) )
				{ continue; }

			if ( edge_array[i].length > 0 )
				{ counter++; }

			parc = edge_array[i].arcs;

			if ( !parc || parc->next )
				{ continue; }

			to_ed = parc->to_ed;
			bal_ed = getTwinEdge ( to_ed );
			parc2 = edge_array[bal_ed].arcs;

			if ( bal_ed == to_ed || !parc2 || parc2->next )
				{ continue; }

			from_ed = i;

			if ( from_ed == to_ed || from_ed == bal_ed )
				{ continue; }

			//linear connection found
			if ( parc->multiplicity <= arcfilter )
			{
				donot1++;
				continue;
			}

			conc_c++;
			bal_fe = getTwinEdge ( from_ed );
			linearUpdateConnection2 ( from_ed, to_ed, 0 );
			allpathUpdateEdge2 ( from_ed, to_ed, 0, last );
			linearUpdateConnection2 ( bal_ed, bal_fe, 1 );
			allpathUpdateEdge2 ( bal_ed, bal_fe, 1, last );
		}

		fprintf ( stderr, "%d edge(s) concatenated in cycle %d.\n", conc_c, round++ );

		if ( arcfilter )
			{ fprintf ( stderr, "%d edge(s) are not linearized because of arc weight is %d.\n", donot1, arcfilter ); }
	}

	fprintf ( stderr, "%d edge(s) in the graph.\n", counter );
}


/*
 * attachPEinfo.c
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
#include "stack.h"
#include "zlib.h"

#define CNBLOCKSIZE 10000
static STACK * isStack;
static int ignorePE1, ignorePE2, ignorePE3, static_flag;
static int onsameCtgPE;
static unsigned long long peSUM;

static boolean staticF;

static int existCounter;

static int calcuIS ( STACK * intStack );

static int cmp_pe ( const void * a, const void * b )
{
	PE_INFO * A, *B;
	A = ( PE_INFO * ) a;
	B = ( PE_INFO * ) b;

	if ( A->rank > B->rank )
	{
		return 1;
	}
	else if ( A->rank == B->rank )
	{
		return 0;
	}
	else
	{
		return -1;
	}
}

/*************************************************
Function:
    loadPEgrads
Description:
    Loads general mapped information of paired-end reads.
Input:
    1. infile:      prefix of output graph file name
Output:
    None.
Return:
    None.
*************************************************/
void loadPEgrads ( char * infile )
{
	FILE * fp;
	char name[256], line[1024];
	int i;
	boolean rankSet = 1;
	sprintf ( name, "%s.peGrads", infile );
	fp = fopen ( name, "r" );

	if ( !fp )
	{
		fprintf ( stderr, "Can not open file %s.\n", name );
		gradsCounter = 0;
		return;
	}

	while ( fgets ( line, sizeof ( line ), fp ) != NULL )
	{
		if ( line[0] == 'g' )
		{
			sscanf ( line + 10, "%d %lld %d", &gradsCounter, &n_solexa, &maxReadLen );
			fprintf ( stderr, "There are %d grad(s), %lld read(s), max read len %d.\n", gradsCounter, n_solexa, maxReadLen );
			break;
		}
	}

	alloc_pe_mem ( gradsCounter );

	for ( i = 0; i < gradsCounter; i++ )
	{
		fgets ( line, sizeof ( line ), fp );
		pes[i].rank = 0;
		sscanf ( line, "%d %lld %d %d", & ( pes[i].insertS ), & ( pes[i].PE_bound ), & ( pes[i].rank ), & ( pes[i].pair_num_cut ) );

		if ( pes[i].rank < 1 )
		{
			rankSet = 0;
		}
	}

	fclose ( fp );

	if ( rankSet )
	{
		qsort ( &pes[0], gradsCounter, sizeof ( PE_INFO ), cmp_pe );
		return;
	}

	int lastRank = 0;

	for ( i = 0; i < gradsCounter; i++ )
	{
		if ( i == 0 )
		{
			pes[i].rank = ++lastRank;
		}
		else if ( pes[i].insertS < 300 )
		{
			pes[i].rank = lastRank;
		}
		else if ( pes[i].insertS < 800 )
		{
			if ( pes[i - 1].insertS < 300 )
			{
				pes[i].rank = ++lastRank;
			}
			else
			{
				pes[i].rank = lastRank;
			}
		}
		else if ( pes[i].insertS < 3000 )
		{
			if ( pes[i - 1].insertS < 800 )
			{
				pes[i].rank = ++lastRank;
			}
			else
			{
				pes[i].rank = lastRank;
			}
		}
		else if ( pes[i].insertS < 7000 )
		{
			if ( pes[i - 1].insertS < 3000 )
			{
				pes[i].rank = ++lastRank;
			}
			else
			{
				pes[i].rank = lastRank;
			}
		}
		else
		{
			if ( pes[i - 1].insertS < 7000 )
			{
				pes[i].rank = ++lastRank;
			}
			else
			{
				pes[i].rank = lastRank;
			}
		}
	}
}

/*************************************************
Function:
    add1Connect
Description:
    Updates the connection between two contigs.
Input:
    1. e1:          left contig
    2. e2:          right contig
    3. gap:         distance between two contigs
    4. weight:      weight of connection
    5. inherit:     inheritance flag of connection
Output:
    None.
Return:
    NULL if failed adding, otherwise the pointer to new or updated connection.
*************************************************/
CONNECT * add1Connect ( unsigned int e1, unsigned int e2, int gap, int weight, boolean inherit )
{
	if ( e1 == e2 || e1 == getTwinCtg ( e2 ) )
	{
		return NULL;
	}

	CONNECT * connect = NULL, *bal_connect = NULL;
	long long sum;

	if ( weight > 255 )
	{
		weight = 255;
	}

	connect = getCntBetween ( e1, e2 );
	bal_connect = getCntBetween ( e2, e1 );

	if ( connect )
	{
		if ( !weight )
		{
			return connect;
		}

		existCounter++;

		if ( !inherit )
		{
			sum = connect->weightNotInherit * connect->gapLen + gap * weight;
			connect->gapLen = sum / ( connect->weightNotInherit + weight );

			if ( connect->weightNotInherit + weight <= 255 )
			{
				connect->weightNotInherit += weight;
			}
			else if ( connect->weightNotInherit < 255 )
			{
				connect->weightNotInherit = 255;
			}
		}
		else
		{
			sum = connect->weight * connect->gapLen + gap * weight;
			connect->gapLen = sum / ( connect->weight + weight );

			if ( !connect->inherit )
			{
				connect->maxSingleWeight = connect->weightNotInherit;
			}

			connect->inherit = 1;
			connect->maxSingleWeight = connect->maxSingleWeight > weight ? connect->maxSingleWeight : weight;
		}

		if ( connect->weight + weight <= 255 )
		{
			connect->weight += weight;
		}
		else if ( connect->weight < 255 )
		{
			connect->weight = 255;
		}
	}
	else
	{
		newCntCounter++;
		connect = allocateCN ( e2, gap );

		if ( cntLookupTable )
		{
			putCnt2LookupTable ( e1, connect );
		}

		connect->weight = weight;

		if ( contig_array[e1].mask || contig_array[e2].mask )
		{
			connect->mask = 1;
		}

		connect->next = contig_array[e1].downwardConnect;
		contig_array[e1].downwardConnect = connect;

		if ( !inherit )
		{
			connect->weightNotInherit = weight;
		}
		else
		{
			connect->weightNotInherit = 0;
			connect->inherit = 1;
			connect->maxSingleWeight = weight;
		}
	}

	return connect;
}

/*************************************************
 Function:
    attach1PE
 Description:
    Checks paired-end relation between two contigs, and updates the connection between
    contigs if the relation is fine.
 Input:
    1. e1:          left contig
    2. pre_pos:     read1's start position on left contig
    3. bal_e2:      right contig in reversed direction
    4. pos:         read2's start position on right contig
    5. insert_size:     insert size of paired-end reads
 Output:
    None.
 Return:
    -1 if left contig and the reversed right contig were the same contig .
    0 if the calculated distance between contigs was abnormal.
    1 if paired-end relation was normal.
    2 if read1 and read2 were aligned to the same contig.
 *************************************************/
int attach1PE ( unsigned int e1, int pre_pos, unsigned int bal_e2, int pos, int insert_size )
{
	int gap, realpeSize;
	unsigned int bal_e1, e2;

	if ( e1 == bal_e2 )
	{
		ignorePE1++;
		return -1;  //orientation wrong
	}

	bal_e1 = getTwinCtg ( e1 );
	e2 = getTwinCtg ( bal_e2 );

	if ( e1 == e2 )
	{
		realpeSize = contig_array[e1].length + overlaplen - pre_pos - pos;

		if ( realpeSize > 0 )
		{
			peSUM += realpeSize;
			onsameCtgPE++;

			if ( ( int ) contig_array[e1].length > insert_size )
			{
				int * item = ( int * ) stackPush ( isStack );
				( *item ) = realpeSize;
			}
		}

		return 2;
	}

	gap = insert_size - overlaplen + pre_pos + pos - contig_array[e1].length - contig_array[e2].length;

	if ( gap < - ( insert_size / 10 ) )
	{
		ignorePE2++;
		return 0;
	}

	if ( gap > insert_size )
	{
		ignorePE3++;
		return 0;
	}

	add1Connect ( e1, e2, gap, 1, 0 );
	add1Connect ( bal_e2, bal_e1, gap, 1, 0 );
	return 1;
}

/*************************************************
 Function:
    connectByPE_grad
 Description:
    Loads alignment information of paired-end reads of specified LIB
    and updates connection between contigs.
 Input:
    1. fp:          alignment information file
    2. peGrad:      LIB number
    3. line:        buffer to store one alignment record
 Output:
    None.
 Return:
    Loaded alignment record number.
 *************************************************/
int connectByPE_grad ( FILE * fp, int peGrad, char * line )
{
	long long pre_readno, readno, minno, maxno;
	int pre_pos, pos, flag, PE, count = 0, Total_PE = 0;
	unsigned int pre_contigno, contigno, newIndex;

	if ( peGrad < 0 || peGrad > gradsCounter )
	{
		fprintf ( stderr, "Specified pe grad is out of bound.\n" );
		return 0;
	}

	maxno = pes[peGrad].PE_bound;

	if ( peGrad == 0 )
	{
		minno = 0;
	}
	else
	{
		minno = pes[peGrad - 1].PE_bound;
	}

	onsameCtgPE = peSUM = 0;
	PE = pes[peGrad].insertS;

	if ( strlen ( line ) )
	{
		sscanf ( line, "%lld %d %d", &pre_readno, &pre_contigno, &pre_pos );

		if ( pre_readno <= minno )
		{
			pre_readno = -1;
		}
	}
	else
	{
		pre_readno = -1;
	}

	ignorePE1 = ignorePE2 = ignorePE3 = 0;
	static_flag = 1;
	isStack = ( STACK * ) createStack ( CNBLOCKSIZE, sizeof ( int ) );

	while ( fgets ( line, lineLen, fp ) != NULL )
	{
		sscanf ( line, "%lld %d %d", &readno, &contigno, &pos );

		if ( readno > maxno )
		{
			break;
		}

		if ( readno <= minno )
		{
			continue;
		}

		newIndex = index_array[contigno];

		if ( isSameAsTwin ( newIndex ) )
		{
			continue;
		}

		if ( PE && ( readno % 2 == 0 ) && ( pre_readno == readno - 1 ) ) // they are a pair of reads
		{
			Total_PE++;
			flag = attach1PE ( pre_contigno, pre_pos, newIndex, pos, PE );

			if ( flag == 1 )
			{
				count++;
			}
		}

		pre_readno = readno;
		pre_contigno = newIndex;
		pre_pos = pos;
	}

	fprintf ( stderr, "For insert size: %d\n", PE );
	fprintf ( stderr, " Total PE links                      %d\n", Total_PE );
	fprintf ( stderr, " Normal PE links on same contig      %d\n", onsameCtgPE );
	fprintf ( stderr, " Abnormal PE links on same contig    %d\n", ignorePE1 );
	fprintf ( stderr, " PE links of minus distance          %d\n", ignorePE2 );
	fprintf ( stderr, " PE links of plus distance           %d\n", ignorePE3 );
	fprintf ( stderr, " Correct PE links                    %d\n", count );
	fprintf ( stderr, " Accumulated connections             %d\n", newCntCounter );
	fprintf ( stderr, "Use contigs longer than %d to estimate insert size: \n", PE );
	fprintf ( stderr, " PE links               %d\n", isStack->item_c );
	calcuIS ( isStack );
	freeStack ( isStack );
	return count;
}

/*************************************************
 Function:
    connectByPE_grad_gz
 Description:
    Loads alignment information of paired-end reads of specified LIB
    and updates connection between contigs.
 Input:
    1. fp:          alignment information file in gz format
    2. peGrad:      LIB number
    3. line:        buffer to store one alignment record
 Output:
    None.
 Return:
    Loaded alignment record number.
 *************************************************/
int connectByPE_grad_gz ( gzFile * fp, int peGrad, char * line )
{
	long long pre_readno, readno, minno, maxno;
	int pre_pos, pos, flag, PE, count = 0, Total_PE = 0;
	unsigned int pre_contigno, contigno, newIndex;

	if ( peGrad < 0 || peGrad > gradsCounter )
	{
		fprintf ( stderr, "Specified pe grad is out of bound.\n" );
		return 0;
	}

	maxno = pes[peGrad].PE_bound;

	if ( peGrad == 0 )
	{
		minno = 0;
	}
	else
	{
		minno = pes[peGrad - 1].PE_bound;
	}

	onsameCtgPE = peSUM = 0;
	PE = pes[peGrad].insertS;

	if ( strlen ( line ) )
	{
		sscanf ( line, "%lld %d %d", &pre_readno, &pre_contigno, &pre_pos );

		if ( pre_readno <= minno )
		{
			pre_readno = -1;
		}
	}
	else
	{
		pre_readno = -1;
	}

	ignorePE1 = ignorePE2 = ignorePE3 = 0;
	static_flag = 1;
	isStack = ( STACK * ) createStack ( CNBLOCKSIZE, sizeof ( int ) );

	while ( gzgets ( fp, line, lineLen ) != NULL )
	{
		sscanf ( line, "%lld %d %d", &readno, &contigno, &pos );

		if ( readno > maxno )
		{
			break;
		}

		if ( readno <= minno )
		{
			continue;
		}

		newIndex = index_array[contigno];

		if ( isSameAsTwin ( newIndex ) )
		{
			continue;
		}

		if ( PE && ( readno % 2 == 0 ) && ( pre_readno == readno - 1 ) ) // they are a pair of reads
		{
			Total_PE++;
			flag = attach1PE ( pre_contigno, pre_pos, newIndex, pos, PE );

			if ( flag == 1 )
			{
				count++;
			}
		}

		pre_readno = readno;
		pre_contigno = newIndex;
		pre_pos = pos;
	}

	fprintf ( stderr, "For insert size: %d\n", PE );
	fprintf ( stderr, " Total PE links                      %d\n", Total_PE );
	fprintf ( stderr, " Normal PE links on same contig      %d\n", onsameCtgPE );
	fprintf ( stderr, " Incorrect oriented PE links         %d\n", ignorePE1 );
	fprintf ( stderr, " PE links of too small insert size   %d\n", ignorePE2 );
	fprintf ( stderr, " PE links of too large insert size   %d\n", ignorePE3 );
	fprintf ( stderr, " Correct PE links                    %d\n", count );
	fprintf ( stderr, " Accumulated connections             %d\n", newCntCounter );
	fprintf ( stderr, "Use contigs longer than %d to estimate insert size: \n", PE );
	fprintf ( stderr, " PE links               %d\n", isStack->item_c );
	calcuIS ( isStack );
	freeStack ( isStack );
	return count;
}

/*************************************************
 Function:
    calcuIS
 Description:
    Calculates insert size of LIB using paired-end reads that both
    reads were aligned to the same contig.
 Input:
    1. intStack:        stack storing the alignment information
 Output:
    None.
 Return:
    Calculated insert size.
 *************************************************/
static int calcuIS ( STACK * intStack )
{
	long long sum = 0;
	int avg = 0;
	int * item;
	int num = intStack->item_c;

	if ( num < 100 )
	{
		fprintf ( stderr, "Too few PE links.\n" );
		return avg;
	}

	stackBackup ( intStack );

	while ( ( item = ( int * ) stackPop ( intStack ) ) != NULL )
	{
		sum += *item;
	}

	stackRecover ( intStack );
	num = intStack->item_c;
	avg = sum / num;
	sum = 0;
	stackBackup ( intStack );

	while ( ( item = ( int * ) stackPop ( intStack ) ) != NULL )
	{
		sum += ( ( long long ) * item - avg ) * ( ( long long ) * item - avg );
	}

	int SD = sqrt ( sum / ( num - 1 ) );

	if ( SD == 0 )
	{
		fprintf ( stderr, " Average insert size    %d\n SD                     %d\n", avg, SD );
		return avg;
	}

	stackRecover ( intStack );
	sum = num = 0;

	while ( ( item = ( int * ) stackPop ( intStack ) ) != NULL )
		if ( abs ( *item - avg ) < 3 * SD )
		{
			sum += *item;
			num++;
		}

	if ( num == 0 ) { avg = 0; }
	else { avg = sum / num; }

	fprintf ( stderr, " Average insert size    %d\n SD                     %d\n", avg, SD );
	return avg;
}

unsigned int getTwinCtg ( unsigned int ctg )
{
	return ctg + contig_array[ctg].bal_edge - 1;
}

boolean isSmallerThanTwin ( unsigned int ctg )
{
	return contig_array[ctg].bal_edge > 1;
}

boolean isLargerThanTwin ( unsigned int ctg )
{
	return contig_array[ctg].bal_edge < 1;
}

boolean isSameAsTwin ( unsigned int ctg )
{
	return contig_array[ctg].bal_edge == 1;
}

/*
 * inc/def.h
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

/* this file provides some datatype definition */
#ifndef _DEF
#define _DEF

#include "def2.h"
#include "types.h"
#include "stack.h"
#include "darray.h"
#include "sam.h"    //support the samfile_t struct

#define EDGE_BIT_SIZE 6
#define word_len 12
#define taskMask 0xf   //the last 7 bits

#define MaxEdgeCov 16000

#define base2int(base)  (char)(((base)&0x06)>>1)    //base ACTG => int 0123
#define int2base(seq)   "ACTG"[seq]                 //int 0123 => base ACTG
#define int2compbase(seq)       "TGAC"[seq]         //int 0123 => base TGAC complement of ACTG
#define int_comp(seq)   (char)(seq^0x02)         //(char)((0x4E>>((seq)<<1))&0x03)

int b_ban;

#ifdef MER127
typedef struct kmer
{
	unsigned long long high1, low1, high2, low2;
} Kmer;
#else
typedef struct kmer
{
	unsigned long long high, low;
} Kmer;
#endif

typedef struct preedge
{
	Kmer from_node;
	Kmer to_node;
	char  *  seq;
	int length;
	unsigned short cvg: 14;
	unsigned bal_edge: 2; //indicate whether it's bal_edge is the previous edge, next edge or itself
} preEDGE;

typedef struct readinterval     //record two paths of bubble
{
	int readid;
	unsigned int edgeid;
	int start;
	struct readinterval * bal_rv;
	struct readinterval * nextOnEdge;   // the downstream in the path
	struct readinterval * prevOnEdge;   // the upstream in the path
	struct readinterval * nextInRead;
	struct readinterval * prevInRead;
} READINTERVAL;

struct arc;
typedef struct edge
{
	unsigned int from_vt;   //from kmer id
	unsigned int to_vt; //to kmer id
	int length;         //edge length
	unsigned short cvg: 14; //coverage
	unsigned short bal_edge: 2; // 2:smaller 0:larger 1:rev-com equal to itself
	unsigned short multi: 14;
	unsigned short deleted : 1;
	unsigned short flag : 1;
	char  *  seq;   //edge content
	READINTERVAL * rv;
	struct arc * arcs;
	long long * markers;    //reads id
} EDGE;

typedef struct edge_sub
{
	unsigned int from_vt;   //from kmer id
	unsigned int to_vt; //to kmer id
	int length;         //edge length
	char  *  seq;   //edge content
} EDGE_SUB;

typedef struct edge_pt
{
	EDGE * edge;
	struct edge_pt * next;
} EDGE_PT;

typedef struct vertex
{
	Kmer kmer;
} VERTEX;
/*
typedef struct connection
{
    unsigned int contigID;
    int gapLen;

    short maxGap;
    unsigned char minGap;
    unsigned char bySmall:1;
    unsigned char weakPoint:1;

    unsigned char weightNotInherit;
    unsigned char weight;
    unsigned char maxSingleWeight;
    unsigned char mask : 1;
    unsigned char used : 1;
    unsigned char weak : 1;
    unsigned char deleted : 1;
    unsigned char prevInScaf : 1;
    unsigned char inherit : 1;
    unsigned char checking : 1;
    unsigned char singleInScaf : 1;
    struct connection *nextInScaf;
    struct connection *next;
    struct connection *nextInLookupTable;
}CONNECT;
*/
typedef struct connection
{
	unsigned int contigID;
	int gapLen;

	unsigned short maxGap;
	unsigned char minGap;
	unsigned char bySmall: 1;
	unsigned char weakPoint: 1;
	unsigned char smallIns: 1;
	unsigned char newIns: 1;

	unsigned char weightNotInherit;
	unsigned char weight;
	unsigned char maxSingleWeight;
	unsigned char mask : 1;
	unsigned char used : 1;
	unsigned char weak : 1;
	unsigned char deleted : 1;
	unsigned char prevInScaf : 1;
	unsigned char inherit : 1;
	unsigned char checking : 1;
	unsigned char singleInScaf : 1;
	struct connection * nextInScaf;
	struct connection * next;
	struct connection * nextInLookupTable;
} CONNECT;

typedef struct prearc
{
	unsigned int to_ed; // the destination edge of prearc
	unsigned int multiplicity;
	struct prearc * next;
} preARC;
/*
typedef struct contig
{
    unsigned int from_vt;
    unsigned int to_vt;
    unsigned int length;
    int to_right;
    unsigned short indexInScaf;
    unsigned char cvg;
    unsigned char bal_edge:2; // 0, 1 or 2
    unsigned char mask : 1;
    unsigned char flag : 1;
    unsigned char multi: 1;
    unsigned char inSubGraph: 1;
    char *seq;
    CONNECT *downwardConnect;
    preARC *arcs;
    STACK *closeReads;
}CONTIG;
*/
typedef struct contig
{
	unsigned int from_vt;       // the first kmer of the contig
	unsigned int to_vt;     // the last kmer of the contig
	unsigned int length;
	unsigned short indexInScaf;     // the index in the scaffold
	unsigned char cvg;
	unsigned char bal_edge: 2; // 0, 1 or 2
	unsigned char mask : 1;
	unsigned char flag : 1;
	unsigned char multi: 1;
	unsigned char inSubGraph: 1;
	unsigned char bubbleInScaff: 1;
	char * seq;
	CONNECT * downwardConnect;      // record the links to other contigs
	preARC * arcs;
	STACK * closeReads;
} CONTIG;

typedef struct read_nearby
{
	int len;
	int dis;  // dis to nearby contig or scaffold's start position
	long long seqStarter;   //sequence start position in dynamic array
} READNEARBY;

typedef struct annotation
{
	unsigned long long readID;
	unsigned int contigID;
	int pos;
} ANNOTATION;

typedef struct parameter
{
	unsigned char threadID;
	void ** hash_table;
	unsigned char * mainSignal;
	unsigned char * selfSignal;
} PARAMETER;

typedef struct lightannot
{
	int contigID;
	int pos;
} LIGHTANNOT;

typedef struct edgepatch
{
	Kmer from_kmer, to_kmer;
	unsigned int length;
	char bal_edge;
} EDGEPATCH;

typedef struct lightctg
{
	unsigned int index;
	int length;
	char  *  seq;
} LIGHTCTG;


typedef struct arc
{
	unsigned int to_ed;
	unsigned int multiplicity;
	struct arc * prev;
	struct arc * next;
	struct arc * bal_arc;
	struct arc * nextInLookupTable;
} ARC;

typedef struct arcexist
{
	Kmer kmer;
	struct arcexist * left;
	struct arcexist * right;
} ARCEXIST;

typedef struct lib_info
{
	int min_ins;
	int max_ins;
	int avg_ins;
	int rd_len_cutoff;  //read length cutoff
	int reverse;
	int asm_flag;
	int map_len;
	int pair_num_cut;
	int rank;
	//indicate which file is next to be read
	int curr_type;
	int curr_index;

	//file handlers to opened files
	FILE * fp1;
	FILE * fp2;
	boolean f1_start;
	boolean f2_start;
	//whether last read is read1 in pair
	int paired;  // 0 -- single; 1 -- read1; 2 -- read2;

	//type1
	char ** a1_fname;
	char ** a2_fname;
	int num_a1_file;
	int num_a2_file;

	//type2
	char ** q1_fname;
	char ** q2_fname;
	int num_q1_file;
	int num_q2_file;

	//type3
	char ** p_fname;
	int num_p_file;  //fasta only

	//type4 &5
	char ** s_a_fname;
	int num_s_a_file;
	char ** s_q_fname;
	int num_s_q_file;

	samfile_t * fp3; //the file handle to read bam file
	char ** b_fname; //the name of the bam file
	int num_b_file; //the number of the bam file
} LIB_INFO;

typedef struct ctg4heap
{
	unsigned int ctgID;
	int dis;
	unsigned char ds_shut4dheap: 1;  // ignore downstream connections
	unsigned char us_shut4dheap: 1;  // ignore upstream connections
	unsigned char ds_shut4uheap: 1;  // ignore downstream connections
	unsigned char us_shut4uheap: 1;  // ignore upstream connections
} CTGinHEAP;

typedef struct ctg4scaf
{
	unsigned int ctgID;
	int start;
	int end;  //position in scaff
	unsigned int cutHead : 8;
	unsigned int cutTail : 7;
	unsigned int scaftig_start : 1;  //is it a scaftig starter
	unsigned int mask : 1;  // is it masked for further operations
	unsigned int gapSeqLen: 15;
	int gapSeqOffset;
} CTGinSCAF;

typedef struct pe_info
{
	int insertS;
	long long PE_bound;
	int rank;
	int pair_num_cut;
} PE_INFO;
#endif

/* this file provides some datatype definition */
#ifndef _DEF
#define _DEF

#include "def2.h"
#include "types.h"
#include "stack.h"
#include "darray.h"

#define EDGE_BIT_SIZE 6
#define word_len 12
#define taskMask 0xf   //the last 7 bits

#define MaxEdgeCov 16000

#define base2int(base)  (char)(((base)&0x06)>>1)
#define int2base(seq)   "ACTG"[seq]
#define int2compbase(seq)       "TGAC"[seq]
#define int_comp(seq)   (char)(seq^0x02) //(char)((0x4E>>((seq)<<1))&0x03)

int b_ban;

typedef unsigned long long Kmer;

typedef struct edon 
{
	Kmer kmer;
	unsigned int ctgLen:1;
	unsigned int twin:1;
	unsigned int pos:30;
	unsigned int ctgID; 
	struct edon *left;
	struct edon *right;
}EDON;

struct node_pt;

typedef struct node
{
	Kmer kmer;
	unsigned char links;
	unsigned char linksB;
	unsigned char cvg;
	unsigned char linear:1;
	unsigned char deleted:1;
	unsigned char mark:1;
	unsigned int to_end;  // the edge no. it belongs to
	struct node *left;
	struct node *right;
}NODE;

typedef struct node_pt
{
	NODE *node;
	Kmer kmer;
	boolean isSmaller;
	struct node_pt *next;
}NODE_PT;

typedef struct preedge
{
	Kmer from_node;
	Kmer to_node;
	char	*seq;
	int length;
	unsigned short cvg;
	unsigned short bal_edge:2;  //indicate whether it's bal_edge is the previous edge, next edge or itself
}preEDGE;

typedef struct readinterval
{
	int readid; 
	unsigned int edgeid;
	int start;
	struct readinterval *bal_rv;
	struct readinterval *nextOnEdge;
	struct readinterval *prevOnEdge;
	struct readinterval *nextInRead;
	struct readinterval *prevInRead;
}READINTERVAL;

struct arc;
typedef struct edge
{
	unsigned int from_vt;
	unsigned int to_vt;
	int length;
	unsigned short cvg:14;
	unsigned short bal_edge:2;
	unsigned short multi:14;
	unsigned short deleted : 1;
	unsigned short flag : 1;
	char	*seq;
	READINTERVAL *rv;
	struct arc *arcs;
	long long *markers;
}EDGE;

typedef struct edge_pt
{
	EDGE *edge;
	struct edge_pt *next;
}EDGE_PT;

typedef struct vertex
{
	Kmer kmer;
}VERTEX;

typedef struct connection
{
	unsigned int contigID;
	int gapLen;
	
	unsigned short maxGap;
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
	int *PE;
}CONNECT;

typedef struct prearc
{
	unsigned int to_ed;
	unsigned int multiplicity;
	struct prearc *next;
}preARC;

typedef struct contig 
{
	unsigned int from_vt;
	unsigned int to_vt;
	unsigned int length;
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

typedef struct read_nearby
{
	int len;
	int dis;  // dis to nearby contig or scaffold's start position
	long long seqStarter;   //sequence start position in dynamic array
}READNEARBY;

typedef struct annotation
{
	unsigned long long readID;
	unsigned int contigID;
	int pos;
}ANNOTATION;

typedef struct parameter
{
	unsigned char threadID;
	void **hash_table;
	unsigned char *mainSignal;
	unsigned char *selfSignal;
}PARAMETER;

typedef struct lightannot
{
	int contigID;
	int pos;
}LIGHTANNOT;

typedef struct edgepatch
{
	Kmer from_kmer,to_kmer;
	unsigned int length;
	char bal_edge;
}EDGEPATCH;

typedef struct lightctg 
{
	unsigned int index;
	int length;
	char	*seq;
}LIGHTCTG;


typedef struct arc
{
	unsigned int to_ed;
	unsigned int multiplicity;
	struct arc *prev;
	struct arc *next;
	struct arc *bal_arc;
	struct arc *nextInLookupTable;
}ARC;

typedef struct arcexist 
{
	Kmer kmer;
	struct arcexist *left;
	struct arcexist *right;
}ARCEXIST;

typedef struct lib_info
{
	int min_ins;
	int max_ins;
	int avg_ins;
	int rd_len_cutoff;
	int reverse;
	int asm_flag;
	int map_len;
	int pair_num_cut;
	int rank;
	//indicate which file is next to be read
	int curr_type;
	int curr_index;
	
	//file handlers to opened files
	FILE *fp1;
	FILE *fp2;
	boolean f1_start;
	boolean f2_start;
	//whether last read is read1 in pair
	int paired;  // 0 -- single; 1 -- read1; 2 -- read2;

//type1
	char **a1_fname;
	char **a2_fname;
	int num_a1_file;
	int num_a2_file;

//type2
	char **q1_fname;
	char **q2_fname;
	int num_q1_file;
	int num_q2_file;

//type3
	char **p_fname;
	int num_p_file;  //fasta only

//type4 &5
	char **s_a_fname;
	int num_s_a_file;
	char **s_q_fname;
	int num_s_q_file;

}LIB_INFO;

typedef struct ctg4heap{
	unsigned int ctgID;
	int dis;
	unsigned char ds_shut4dheap:1;   // ignore downstream connections
	unsigned char us_shut4dheap:1;   // ignore upstream connections
	unsigned char ds_shut4uheap:1;   // ignore downstream connections
	unsigned char us_shut4uheap:1;   // ignore upstream connections
}CTGinHEAP;

typedef struct ctg4scaf{
	unsigned int ctgID;
	int start;
	int end;  //position in scaff
	unsigned int cutHead : 8; //
	unsigned int cutTail : 7; //
	unsigned int scaftig_start : 1;  //is it a scaftig starter
	unsigned int mask : 1;  // is it masked for further operations 
	unsigned int gapSeqLen:15;
	int gapSeqOffset;
}CTGinSCAF;

typedef struct pe_info{
	int insertS;
	long long PE_bound;
	int rank;
	int pair_num_cut;
}PE_INFO;
#endif

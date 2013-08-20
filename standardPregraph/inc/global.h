/*
 * inc/global.h
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

int visual = 0;         // 1 for output some files , which are useful for visual
int * contig_index_array = NULL;
int scaffNum = 0;
int gapNum = 1;
boolean fill = 0; // 1 for output some files ,which are useful for the software "kgf"
int overlaplen = 23;//k-mer Size
int inGraph;        //for checking whether -g  is set, (graph prefix)
long long n_ban;    //not used
long long n_solexa = 0; //reads number
long long prevNum = 0;  //not used
int ins_size_var = 20;      // SD of insert-size
PE_INFO * pes = NULL;     //record the pe info in lib file
MEM_MANAGER * rv_mem_manager = NULL;
MEM_MANAGER * cn_mem_manager = NULL;
MEM_MANAGER * arc_mem_manager = NULL;
unsigned int num_vt = 0;                // num of the end-kmer
unsigned long long new_num_vt = 0;      // the new num of the end-kmer after adding the new end-kmer
unsigned int ** found_routes = NULL;
unsigned int * so_far = NULL;       // recorf the path of contig while filling gap
int max_n_routes = 10;
int num_trace;
Kmer WORDFILTER;        //mask code for extracting Kmer info from raw data (two unsigned long long int)
unsigned int num_ed = 0;  //number of edges
unsigned int num_ctg = 0;               // num of contig
unsigned int num_ed_limit;              // the count of edge
unsigned int extraEdgeNum;          // the new count of edge after adding the new edge
EDGE * edge_array = NULL;           // used to record all the info of edge
VERTEX * vt_array = NULL;               // used to record the sequence info of the end-kmer
unsigned int * index_array = NULL;      // used to translate the old contig index to the new contig index
CONTIG * contig_array = NULL;       // used to record all the info of contig
int lineLen;
int len_bar = 100;
int weakPE = 3;     // the minimun weight requirement for the connection
int fillGap = 0;        // 1 for fill the gap after scaffold asm
boolean globalFlag;
long long arcCounter;       // record the num of the arc
MEM_MANAGER * prearc_mem_manager = NULL;
MEM_MANAGER ** preArc_mem_managers = NULL;
int maxReadLen = 0;       //max length will be used for each LIB, soapdenovo read LIBs one by one , for each set a maxReadLen
int maxReadLen4all = 0;   //max length will be used for all reads
int minReadLen = 0;     // min length will be used for all readss
int maxNameLen = 0;       //max length for the name of reads or sequences
ARC ** arcLookupTable = NULL;
long long * markersArray = NULL;
boolean deLowKmer = 0;  //remove the kmers which coverage are not bigger than deLowKmer
boolean deLowEdge = 1;  //remove the edges which coverage are not bigger than deLowEdge
long long newCntCounter;        // record the number of the new connection in one insert-size
long long discardCntCounter;
boolean repsTie = 0;            //sovle tiny repeat or not
CONNECT ** cntLookupTable = NULL;
int num_libs = 0;                   //number of LIBs in read config file
LIB_INFO * lib_array = NULL;        //store LIB's info into lib_array
int libNo = 0;                  // the current number of lib
long long readNumBack;
int gradsCounter;                  //pair number in lib file
unsigned int ctg_short = 0;     //shortest contig for scaffolding
int thrd_num = 8;                  //thread number
int cvgAvg = 0; // the average coverage of contigs
KmerSet ** KmerSets = NULL;      //KmerSet [i] for  thread i
KmerSet ** KmerSetsPatch = NULL; //KmerSet for (k+1) mer
DARRAY * readSeqInGap = NULL;
DARRAY * gapSeqDarray = NULL;
DARRAY ** darrayBuf;
boolean orig2new;   // 1 for re-arrange the contig index using the length
int maxSteps;
boolean maskRep = 1;        // 1 for masking repeat for scaffold asm , 0 for un-masking repeat.
int GLDiff = 50;
int initKmerSetSize = 0;   // init_size = (ubyte8) ((double) initKmerSetSize * 1024.0f * 1024.0f * 1024.0f / (double) thrd_num / 24);
long known_genome_size = 0;
int smallKmer = 0;  // the kmer of the step "Map"
int deltaKmer = 0;  // for map, K-k

double cvg_low = 0.1;
double cvg_high = 2;
double len_times = 0;

float ins_var_idx = 1.5;
int Insert_size = 0;        // the current insert-size
int score_mask = 1;
int COMPATIBLE_MODE = 0;     // 1 for the gz file ; 0 for the normal file
float cvg4SNP = 0.6;

MEM_MANAGER * edgeid_mem_manager = NULL;
unsigned int num_vtnew = 0; //new vertex num
unsigned int kmer_cnew = 0; //new kmer num
const int step = 1;     //step for multi kmer
//int nowstep = 1;
int nowstep2 = 1;

unsigned int * edge_id = NULL;  //edge id array
VERTEX * vt_arraynew = NULL;    //vertex array for k+1mer

KmerSet2 * KmerSetsNew = NULL; //kmer set for k+1mer
char libfilename[256];
boolean parse = 0;
unsigned int num_ed_temp = 0;   //  record the count of the edge

int arcfilter = 0;  //arc filter thrd
boolean outputContig = 0;

long long pinCounter;       //the count of the merged bubble
int clean;  //merge clean bubble

unsigned int num_kmer_limit;

int gLineLen = 5000;
char * gStr = NULL;

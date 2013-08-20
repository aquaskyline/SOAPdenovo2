/*
 * inc/extvab.h
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

/*** global variables ****/
extern boolean fill;
extern int overlaplen;
extern int inGraph;
extern long long n_ban;
extern Kmer WORDFILTER;
extern boolean globalFlag;
extern int thrd_num;
extern int gLineLen;
extern char * gStr;
/**** reads info *****/
extern long long n_solexa;
extern long long prevNum;
extern int ins_size_var;
extern PE_INFO * pes;
extern int maxReadLen;
extern int maxReadLen4all;
extern int minReadLen;
extern int maxNameLen;
extern int num_libs;
extern LIB_INFO * lib_array;
extern int libNo;
extern long long readNumBack;
extern int gradsCounter;
/*** used for pregraph *****/
extern MEM_MANAGER * prearc_mem_manager; //also used in scaffolding
extern MEM_MANAGER ** preArc_mem_managers;
extern boolean deLowKmer;
extern boolean deLowEdge;
extern KmerSet ** KmerSets; // also used in mapping
extern KmerSet ** KmerSetsPatch;
/**** used for contiging ****/
extern boolean repsTie;
extern long long arcCounter;
extern unsigned int num_ed;
extern unsigned int num_ed_limit;
extern unsigned int extraEdgeNum;
extern EDGE * edge_array;
extern VERTEX * vt_array;
extern MEM_MANAGER * rv_mem_manager;
extern MEM_MANAGER * arc_mem_manager;
extern unsigned int num_vt;
extern unsigned long long new_num_vt;
extern int len_bar;
extern ARC ** arcLookupTable;
extern long long * markersArray;
/***** used for scaffolding *****/
extern MEM_MANAGER * cn_mem_manager;
extern unsigned int num_ctg;
extern unsigned int * index_array;
extern CONTIG * contig_array;
extern int lineLen;
extern int weakPE;
extern long long newCntCounter;
extern CONNECT ** cntLookupTable;
extern unsigned int ctg_short;
extern int cvgAvg;
extern boolean orig2new;
extern int Insert_size;
extern int score_mask;
extern long long discardCntCounter;
extern int COMPATIBLE_MODE;
extern float cvg4SNP;
/**** used for gapFilling ****/
extern DARRAY * readSeqInGap;
extern DARRAY * gapSeqDarray;
extern DARRAY ** darrayBuf;
extern int fillGap;
/**** used for searchPath *****/
extern int maxSteps;
extern int num_trace;
extern unsigned int ** found_routes;
extern unsigned int * so_far;
extern int max_n_routes;
extern boolean maskRep;
extern int GLDiff;
extern int initKmerSetSize;
extern long known_genome_size;
extern int smallKmer;
extern int deltaKmer;
extern int gapNum;
extern int scaffNum;
extern int * contig_index_array;

extern double cvg_low;
extern double cvg_high;
extern double len_times;
extern float ins_var_idx;
extern int visual;

extern unsigned int num_vtnew;
extern unsigned int kmer_cnew;
extern const int step;
extern unsigned int * edge_id;
extern VERTEX * vt_arraynew;
extern KmerSet2 * KmerSetsNew;
extern MEM_MANAGER * edgeid_mem_manager;
extern char libfilename[256];
extern boolean parse;
extern int nowstep2;
extern unsigned int num_ed_temp;

extern int arcfilter;
extern boolean outputContig;

extern int clean;
extern long long pinCounter;
extern unsigned int num_kmer_limit;


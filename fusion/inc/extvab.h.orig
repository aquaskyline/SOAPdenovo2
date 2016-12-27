/*************************************************************************** 
 * Title:          extvab.h 
 * Author:         Hongmei Zhu 
 * Created:        Jun. 2007 
 * Last modified:  May. 2009 
 * 
 * All Rights Reserved 
 * See file LICENSE for details. 
 ***************************************************************************/ 
/*** global variables ****/
extern int overlaplen;
extern int inGraph;
extern long long n_ban; 
extern Kmer WORDFILTER;
extern boolean globalFlag;
extern int thrd_num;

extern int verbosity;
extern char verboseStr[verboseBufSize];

/**** reads info *****/
extern long long n_solexa; 
extern long long prevNum; 
extern int ins_size_var;
extern PE_INFO *pes;
extern int maxReadLen;
extern int maxReadLen4all;
extern int minReadLen;
extern int maxNameLen;
extern int num_libs;
extern LIB_INFO *lib_array;
extern int libNo;
extern long long readNumBack;
extern int gradsCounter;
/*** used for pregraph *****/
extern MEM_MANAGER *prearc_mem_manager;  //also used in scaffolding
extern MEM_MANAGER **preArc_mem_managers;
extern boolean deLowKmer;
extern boolean deLowEdge;
extern KmerSet **KmerSets;  // also used in mapping
extern KmerSet **KmerSetsPatch;

extern spcKmerSet *spcSet;

/**** used for contiging ****/
extern boolean repsTie;
extern long long arcCounter;
extern unsigned int num_ed;
extern unsigned int num_ed_limit;
extern unsigned int extraEdgeNum;
extern EDGE *edge_array;
extern VERTEX *vt_array;
extern MEM_MANAGER *rv_mem_manager;
extern MEM_MANAGER *arc_mem_manager;
extern unsigned int num_vt;
extern int len_bar;
extern ARC **arcLookupTable;
extern long long *markersArray;
/***** used for scaffolding *****/
extern MEM_MANAGER *cn_mem_manager;
extern unsigned int num_ctg;
extern unsigned int *index_array;
extern CONTIG *contig_array;
extern int lineLen;
extern int weakPE;
extern long long newCntCounter;
extern CONNECT **cntLookupTable;
extern unsigned int ctg_short;
extern int cvgAvg;
extern boolean orig2new;
/**** used for gapFilling ****/
extern DARRAY *readSeqInGap;
extern DARRAY *gapSeqDarray;
extern DARRAY **darrayBuf;
extern int fillGap;
/**** used for searchPath *****/
extern int maxSteps;
extern int num_trace;
extern unsigned int**found_routes;
extern unsigned int*so_far;
extern int max_n_routes;
extern boolean maskRep;
extern int GLDiff;
extern int initKmerSetSize;
extern char *shortrdsfile;
extern char *graphfile;
extern double OverlapPercent ;
extern double ConflPercent ;
extern double close_threshold;
extern int bund_threshold;
extern char *ctg_file;
//extern boolean large_kmer;

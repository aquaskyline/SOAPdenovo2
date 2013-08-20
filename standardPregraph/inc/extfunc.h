/*
 * inc/extfunc.h
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

#include "check.h"
#include "extfunc2.h"

extern void initAIO ( struct aiocb * aio, char * buf, int fd, int size );
extern int AIORead ( struct aiocb * mycb, int * offset, char * buf, char * cach, int * rt, int curr_type );
extern boolean check_file ( char * name ); //add 2012.7.6
extern boolean checkFiles4Scaff ( char * infile );

extern boolean openNextFile ( int * libNo, boolean pairs, unsigned char asm_ctg );
extern int nextValidIndex ( int libNo, boolean pair, unsigned char asm_ctg );
extern void openFileInLib ( int libNo );
extern void closeFp1InLab ( int libNo );
extern void closeFp2InLab ( int libNo );
extern boolean readseqInLib ( char * src_seq, char * src_name, int * len_seq, char * buf, int * start, int offset, int i );

extern void readseq1by1 ( char * src_seq, char * src_name, int * len_seq, FILE * fp, long long num_seq );
extern void readseqPbyP ( char * src_seq, char * src_name, int * insertS, int * len_seq, FILE * fp, long long num_seq );
extern long long readseqpar ( int * max_len, int * min_leg, int * max_name_len, FILE * fp );
extern void free_edge_list ( EDGE_PT * el );
extern void reverseComplementSeq ( char * seq, int len, char * bal_seq );
extern void free_edge_array ( EDGE * ed_array, int ed_num );
extern void free_lightctg_array ( LIGHTCTG * ed_array, int ed_num );
extern char getCharInTightString ( char * tightSeq, int pos );
extern void writeChar2tightSting ( char nt, char * tightSeq, int pos );
extern void short_reads_sum();
extern void read_one_sequence ( FILE * fp, long long * T, char ** X );
extern void output_edges ( preEDGE * ed_array, int ed_num, char * outfile );
extern void loadVertex ( char * graphfile );
extern void loadEdge ( char * graphfile );
extern boolean loadPath ( char * graphfile );
extern READINTERVAL * allocateRV ( int readid, int edgeid );
extern void createRVmemo();
extern void dismissRV ( READINTERVAL * rv );
extern void destroyReadIntervMem();
extern void destroyConnectMem();
extern void u2uConcatenate();
extern void output_contig ( EDGE * ed_array, unsigned int ed_num, char * outfile, int cut_len );
extern void printTightString ( char * tightSeq, int len );
extern int roughUniqueness ( unsigned int edgeno, char ignore_cvg, char * ignored );
extern void outputReadPos ( char * graphfile, int min_len );
extern void testSearch();
extern void allpathConcatenate();
extern void output_updated_edges ( char * outfile );
extern void output_updated_vertex ( char * outfile );
extern void loadUpdatedEdges ( char * graphfile );
extern void loadUpdatedVertex ( char * graphfile );
extern void connectByPE ( char * infile );
extern void output_cntGVZ ( char * outfile );
extern void output_graph ( char * outfile );
extern void testLinearC2C();
extern void output_contig_graph ( char * outfile );
extern void scaffolding ( unsigned int cut_len, char * outfile );
extern int cmp_int ( const void * a, const void * b );
extern CONNECT * allocateCN ( unsigned int contigId, int gap );
extern int recoverRep();
extern void loadPEgrads ( char * infile );
extern int putInsertS ( long long readid, int size, int * currGrads );
extern int getInsertS ( long long readid, int * readlen );
extern int connectByPE_grad ( FILE * fp, int peGrad, char * line );
extern int connectByPE_grad_gz ( gzFile * fp, int peGrad, char * line );
extern void PEgradsScaf ( char * infile );
extern void reorderAnnotation ( char * infile, char * outfile );
extern void output_1edge ( preEDGE * edge, gzFile * fp );
extern void prlRead2edge ( char * libfile, char * outfile );
extern void annotFileTrans ( char * infile, char * outfile );
extern void prlLoadPath ( char * graphfile );
extern void misCheck ( char * infile, char * outfile );
extern int uniqueLenSearch ( unsigned int * len_array, unsigned int * flag_array, int num, unsigned int target );
extern int cmp_vertex ( const void * a, const void * b );
extern void linkContig2Vts();
extern int connectByPE_gradPatch ( FILE * fp1, FILE * fp2, int peGrad, char * line1, char * line2 );
extern void scaftiging ( char * graphfile, int len_cut );
extern void gapFilling ( char * graphfile, int cut_len );
extern ARC * getArcBetween ( unsigned int from_ed, unsigned int to_ed );
extern void bubblePinch ( double simiCutoff, char * outfile, int M, boolean isIter, boolean last );
extern void linearConcatenate ( boolean isIter, boolean last );
extern unsigned char setArcMulti ( unsigned int from_ed, unsigned int to_ed, unsigned char value );
extern ARC * allocateArc ( unsigned int edgeid );
extern void cutTipsInGraph ( int cutLen, boolean strict, boolean last );
extern ARC * deleteArc ( ARC * arc_list, ARC * arc );
extern void compactEdgeArray();
extern void dismissArc ( ARC * arc );
extern void createArcMemo();
extern ARC * getArcBetween ( unsigned int from_ed, unsigned int to_ed );
extern ARC * allocateArc ( unsigned int edgeid );
extern void writeChar2tightString ( char nt, char * tightSeq, int pos );
extern void output_heavyArcs ( char * outfile );
extern preARC * allocatePreArc ( unsigned int edgeid );
extern void destroyPreArcMem();
extern void traceAlongArc ( unsigned int destE, unsigned int currE, int max_steps, int min, int max, int index, int len, int * num_route );
extern void freeContig_array();
extern void output_scafSeq ( char * graphfile, int len_cut );
extern void putArcInHash ( unsigned int from_ed, unsigned int to_ed );
extern boolean DoesArcExist ( unsigned int from_ed, unsigned int to_ed );
extern void recordArcInHash();
extern void destroyArcHash();
extern void removeWeakEdges ( int lenCutoff, unsigned int multiCutoff );
extern void createArcLookupTable();
extern void deleteArcLookupTable();
extern void putArc2LookupTable ( unsigned int from_ed, ARC * arc );
extern void removeArcInLookupTable ( unsigned int from_ed, unsigned int to_ed );
extern ARC * arcCount ( unsigned int edgeid, unsigned int * num );
extern void mapFileTrans ( char * infile );
extern void solveReps();
extern void removeDeadArcs();
extern void destroyArcMem();
extern void getCntsInFile ( char * infile );
extern void scafByCntInfo ( char * infile );
extern CONNECT * add1Connect ( unsigned int e1, unsigned int e2, int gap, int weight, boolean inherit );
extern void getScaff ( char * infile );
extern void traceAlongMaskedCnt ( unsigned int destE, unsigned int currE, int max_steps, int min, int max, int index, int len, int * num_route );
extern void createPreArcMemManager();
extern boolean loadPathBin ( char * graphfile );
extern void recordArcsInLookupTable();
extern FILE * multiFileRead1seq ( char * src_seq, char * src_name, int * len_seq, FILE * fp, FILE * freads );
extern void multiFileSeqpar ( FILE * fp );
extern long long multiFileParse ( int * max_leg, int * min_leg, int * max_name_leg, FILE * fp );
extern CONNECT * getCntBetween ( unsigned int from_ed, unsigned int to_ed );
extern void createCntMemManager();
extern void destroyConnectMem();
extern void createCntLookupTable();
extern void deleteCntLookupTable();
extern void putCnt2LookupTable ( unsigned int from_c, CONNECT * cnt );
extern void prlRead2Ctg ( char * seqfile, char * outfile );
extern boolean prlContig2nodes ( char * grapfile, int len_cut );
extern void scan_libInfo ( char * libfile );
extern void free_libs();

extern boolean read1seqInLibBam ( char * src_seq, char * src_name, int * len_seq, int * libNo, boolean pair, unsigned char asm_ctg, int * type );
extern boolean read1seqInLib ( char * src_seq, char * src_name, int * len_seq,
                               int * libNo, boolean pair, unsigned char asm_ctg , int * type );
extern void save4laterSolve();
extern void solveRepsAfter();
extern void free_pe_mem();
extern void alloc_pe_mem ( int gradsCounter );
extern void prlDestroyPreArcMem();
extern preARC * prlAllocatePreArc ( unsigned int edgeid, MEM_MANAGER * manager );
extern boolean prlRead2HashTable ( char * libfile, char * outfile );
extern void free_allSets();
extern void removeSingleTips();
extern void removeMinorTips();
extern void kmer2edges ( char * outfile );
extern void output_vertex ( char * outfile );
extern boolean prlRead2HashTable ( char * libfile, char * outfile );
extern void Links2Scaf ( char * infile );
extern void PE2Links ( char * infile );
extern unsigned int getTwinCtg ( unsigned int ctg );
extern void basicContigInfo ( char * infile );
extern boolean isSmallerThanTwin ( unsigned int ctg );
extern boolean isLargerThanTwin ( unsigned int ctg );
extern boolean isSameAsTwin ( unsigned int ctg );
extern boolean loadMarkerBin ( char * graphfile );
extern void readsCloseGap ( char * graphfile );
extern void prlReadsCloseGap ( char * graphfile );
extern void locateReadOnScaf ( char * graphfile );
/*********** Kmer related *************/
extern Kmer createFilter ( int overlaplen );
extern void printKmerSeq ( FILE * fp, Kmer kmer );
//extern U256b Kmer2int256(Kmer seq);
extern boolean KmerLarger ( Kmer kmer1, Kmer kmer2 );
extern boolean KmerSmaller ( Kmer kmer1, Kmer kmer2 );
extern boolean KmerEqual ( Kmer kmer1, Kmer kmer2 );
extern Kmer KmerAnd ( Kmer kmer1, Kmer kmer2 );
extern Kmer KmerLeftBitMoveBy2 ( Kmer word );
extern Kmer KmerRightBitMoveBy2 ( Kmer word );
extern Kmer KmerPlus ( Kmer prev, char ch );
extern Kmer nextKmer ( Kmer prev, char ch );
extern Kmer prevKmer ( Kmer next, char ch );
extern char firstCharInKmer ( Kmer kmer );
extern Kmer KmerRightBitMove ( Kmer word, int dis );
extern Kmer reverseComplement ( Kmer word, int overlap );
extern ubyte8 hash_kmer ( Kmer kmer );
extern int kmer2vt ( Kmer kmer );
extern void print_kmer ( FILE * fp, Kmer kmer, char c );
extern int bisearch ( VERTEX * vts, int num, Kmer target );
extern void printKmerSeq ( FILE * fp, Kmer kmer );
extern char lastCharInKmer ( Kmer kmer );
int localGraph ( READNEARBY * rdArray, int num, CTGinSCAF * ctg1, CTGinSCAF * ctg2,
                 int origOverlap, Kmer * kmerCtg1, Kmer * kmerCtg2,
                 int overlap, DARRAY * gapSeqArray, char * seqCtg1, char * seqCtg2, char * seqGap );
extern unsigned int getTwinEdge ( unsigned int edgeno );
extern boolean EdSmallerThanTwin ( unsigned int edgeno );
extern boolean EdLargerThanTwin ( unsigned int edgeno );
extern boolean EdSameAsTwin ( unsigned int edgeno );
extern void removeLowCovEdges ( int lenCutoff, unsigned short covCutoff, boolean last );
extern int getMaxLongReadLen ( int num_libs );
extern void prlLongRead2Ctg ( char * libfile, char * outfile );
extern void outputTightStr ( FILE * fp, char * tightStr, int start, int length, int outputlen, int revS, int * col );
extern void crc32c_Init();

extern int validArcCount ( preARC * arc, int cutoff );
extern unsigned int maxArcWeight ( preARC * arc );
extern __uint128_t Kmer2int128 ( Kmer seq );
extern void printSeq ( FILE * fo, char * seq, int len );


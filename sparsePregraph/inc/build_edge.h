/*
 * inc/sparse_kmer.h
 *
 * Copyright (c) 2008-2016 Ruibang Luo <aquaskyline.com>.
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

#ifndef _BUILD_EDGE_H

#define _BUILD_EDGE_H

void removeMinorTips ( struct hashtable2 *ht, int K_size, int cut_len_tip, int &tip_c );
void kmer2edges ( hashtable2 *ht, int K_size, char *outfile );
void convert ( char *sparse_edge_file, int K_size, char *output_prefix );
void RemovingWeakNodesAndEdges2 ( hashtable2 *ht, int K_size, int NodeCovTh, int EdgeCovTh, size_t *bucket_cnt, size_t *edge_cnt );

struct stacked_node2
{
  struct bucket2 *node;
  bool is_left; // change it to a byte later
  struct edge_node *edge;
  struct stacked_node2 *next;
};

typedef struct preedge2
{
  struct stacked_node2 *from_node;
  struct stacked_node2 *to_node;
  string *full_edge;
  unsigned short cvg;
  unsigned short bal_edge: 2;

} preEDGE2;


// below is static methods
//remove minor tips ...
//void removeMinorTips(struct hashtable2 *ht,int K_size,int cut_len_tip,int &tip_c);
static void mask1in1out ( hashtable2 *ht );
static int clipTipFromNode ( hashtable2 *ht, int K_size, bucket2 *node, int cut_len_tip );
static int count_left_edge_num ( bucket2 *bkt );
static int count_right_edge_num ( bucket2 *bkt );

static void dislink ( hashtable2 *ht, int K_size, stacked_node2 *from_node );
static bucket2 *lastKmer ( hashtable2 *ht, int K_size, bucket2 *node, edge_node *edge, int is_left, int &smaller );
//static bucket2* search_kmer(hashtable2 *ht,uint64_t* t_kmer, int Kmer_arr_sz); old
static bucket2 *search_kmer ( hashtable2 *ht, kmer_t2 *t_kmer );

static void removeEdge ( bucket2 *node, edge_node *edge, int is_left );
static void stat_edge_num ( hashtable2 *ht );
static void stat_edge_cvg_len ( hashtable2 *ht );
static bool isSmaller2 ( uint64_t *kmer, int K_size );

//kmer2edges ....
//void kmer2edges(hashtable2* ht,int K_size,char *outfile);
static void make_edge ( hashtable2 *ht, int K_size, FILE *fp );
static int startEdgeFromNode ( hashtable2 *ht, int K_size, bucket2 *node, FILE *fp );
static void stringBeads ( hashtable2 *ht, int K_size, list<stacked_node2 *> &stack, stacked_node2 *from_node, edge_node *from_edge, int *node_c );

static void process_1stack ( hashtable2 *ht, int K_size, list<stacked_node2 *> &stack, FILE *fp, vector<preEDGE2> &loops_edges );

//static void  get_kmer(const char * seq,int len, int K_size,int pos,uint64_t *kmer,int arr_sz );
static void output_1edge ( preEDGE2 *long_edge, int K_size, FILE *fp );
static string stack2string ( hashtable2 *ht, int K_size, list<stacked_node2 *> &stack );
static bool check_palindrome ( string &str );
static string revCompSeq ( const string &str );

//convert the edge fomat ...
//void convert(char * sparse_edge_file,int K_size, char * output_prefix);
static void convert_kmer ( uint64_t *sparse_kmer, int K_size, int arr_sz );
static  uint64_t *fastReverseComp ( uint64_t *seq_arr, int seq_size, int arr_sz );

#endif









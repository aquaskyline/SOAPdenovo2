/*
 * inc/newhash.h
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

#ifndef __NEW_HASH_RJ
#define __NEW_HASH_RJ

#ifndef K_LOAD_FACTOR
#define K_LOAD_FACTOR 0.75
#endif

#define MAX_KMER_COV 63
#define EDGE_BIT_SIZE 6 //6 bit for each edge
#define EDGE_XOR_MASK 0x3FU
#define LINKS_BITS 0x00FFFFFFU

#define get_kmer_seq(mer) ((mer).seq)
#define set_kmer_seq(mer, val) ((mer).seq = val)

#define get_kmer_left_cov(mer, idx) (((mer).l_links>>((idx)*EDGE_BIT_SIZE))&EDGE_XOR_MASK)
#define set_kmer_left_cov(mer, idx, val) ((mer).l_links = ((mer).l_links&(~(EDGE_XOR_MASK<<((idx)*EDGE_BIT_SIZE)))) | (((val)&EDGE_XOR_MASK)<<((idx)*EDGE_BIT_SIZE)) )
#define get_kmer_left_covs(mer) (get_kmer_left_cov(mer, 0) + get_kmer_left_cov(mer, 1) + get_kmer_left_cov(mer, 2) + get_kmer_left_cov(mer, 3))

#define get_kmer_right_cov(mer, idx) (((mer).r_links>>((idx)*EDGE_BIT_SIZE))&EDGE_XOR_MASK)
#define set_kmer_right_cov(mer, idx, val) ((mer).r_links = ((mer).r_links&(~(EDGE_XOR_MASK<<((idx)*EDGE_BIT_SIZE)))) | (((val)&EDGE_XOR_MASK)<<((idx)*EDGE_BIT_SIZE)) )
#define get_kmer_right_covs(mer) (get_kmer_right_cov(mer, 0) + get_kmer_right_cov(mer, 1) + get_kmer_right_cov(mer, 2) + get_kmer_right_cov(mer, 3))


#define is_kmer_entity_null(flags, idx)    ((flags)[(idx)>>4]>>(((idx)&0x0f)<<1)&0x01)
#define is_kmer_entity_del(flags, idx)     ((flags)[(idx)>>4]>>(((idx)&0x0f)<<1)&0x02)
#define set_kmer_entity_null(flags, idx)   ((flags)[(idx)>>4] |= (0x01u<<(((idx)&0x0f)<<1)))
#define set_kmer_entity_del(flags, idx)    ((flags)[(idx)>>4] |= (0x02u<<(((idx)&0x0f)<<1)))
#define clear_kmer_entity_null(flags, idx) ((flags)[(idx)>>4] &= ~(0x01u<<(((idx)&0x0f)<<1)))
#define clear_kmer_entity_del(flags, idx)  ((flags)[(idx)>>4] &= ~(0x02u<<(((idx)&0x0f)<<1)))
#define exists_kmer_entity(flags, idx)     (!((flags)[(idx)>>4]>>(((idx)&0x0f)<<1)&0x03))

#ifdef MER127
typedef __uint128_t u128b;

typedef struct u256b
{
	u128b low;
	u128b high;
} U256b;
#else
#endif

typedef struct kmer_st
{
	Kmer seq;       //kmer set
	ubyte4 l_links;                     // sever as edgeID since make_edge
	ubyte4 r_links: 4 * EDGE_BIT_SIZE;
	ubyte4 linear: 1;
	ubyte4 deleted: 1;
	ubyte4 checked: 1;
	ubyte4 single: 1;
	ubyte4 twin: 2;
	ubyte4 inEdge: 2;
} kmer_t;

typedef struct kmerSet_st
{
	kmer_t * array; //kmer set
	ubyte4 * flags; //mark the element pos that exist in array
	ubyte8 size;
	ubyte8 count;
	ubyte8 max;
	float load_factor;
	ubyte8 iter_ptr;    //iter in set
} KmerSet;

typedef struct kmer_pt
{
	kmer_t * node;
	Kmer kmer;
	boolean isSmaller;
	struct kmer_pt * next;
} KMER_PT;

extern KmerSet * init_kmerset ( ubyte8 init_size, float load_factor );
extern int search_kmerset ( KmerSet * set, Kmer seq, kmer_t ** rs );
extern int put_kmerset ( KmerSet * set, Kmer seq, ubyte left, ubyte right, kmer_t ** kmer_p );
extern byte8 count_kmerset ( KmerSet * set );
extern void free_Sets ( KmerSet ** KmerSets, int num );
extern void free_kmerset ( KmerSet * set );
extern void dislink2nextUncertain ( kmer_t * node, char ch, boolean smaller );
extern void dislink2prevUncertain ( kmer_t * node, char ch, boolean smaller );

extern int count_branch2prev ( kmer_t * node );
extern int count_branch2next ( kmer_t * node );
extern char firstCharInKmer ( Kmer kmer );

#endif

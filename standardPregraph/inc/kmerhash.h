/*
 * inc/kmerhash.h
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

#ifndef __KMER_HASH_RJ
#define __KMER_HASH_RJ

#define set_kmer_id(mer, val) ((mer).edgeId = val)

typedef struct edgeID
{
	unsigned int edge;      //edge id
	char flag;  //00: big to, 01: big from, 10: small to, 11: small from
	struct edgeID * next;
} EDGEID;

typedef struct kmer_st2
{
	Kmer seq;   //kmer sequence
	ubyte4 l_links; //left out degree
	ubyte4 r_links; //right out degree
	int count;  //edge number
	struct edgeID * edgeId;
} kmer_t2;

typedef struct kmerSet_st2
{
	kmer_t2 * array; //kmer set
	ubyte4 * flags; //mark the element pos that exist in array
	ubyte8 size;
	ubyte8 count;
	ubyte8 max;
	double load_factor;
	ubyte8 iter_ptr;
} KmerSet2;

extern KmerSet2 * init_kmerset2 ( ubyte8 init_size, float load_factor );
extern int search_kmerset2 ( KmerSet2 * set, Kmer seq, kmer_t2 ** rs );
extern int put_kmerset2 ( KmerSet2 * set, Kmer seq, int id, char flag, kmer_t2 ** kmer_p );
extern byte8 count_kmerset2 ( KmerSet2 * set );
extern void free_Sets2 ( KmerSet2 ** KmerSets, int num );
extern void free_kmerset2 ( KmerSet2 * set );
extern void update_kmer2 ( kmer_t2 * mer, int id, char flag );
extern void set_new_kmer2 ( kmer_t2 * mer, Kmer seq, int id, char flag );

#endif


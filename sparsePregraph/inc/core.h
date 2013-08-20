/*
 * inc/core.h
 *
 * This file is part of SOAPdenovo.
 *
 * Part of this file is refered and modified from SparseAssembler
 * (See <http://sourceforge.net/projects/sparseassembler/>).
 *
 */

#ifndef _CORE_H
#define _CORE_H

#include "stdinc.h"


#ifdef _63MER_
struct kmer_t2 //use union later
{
	uint64_t kmer[2];
};
#endif


#ifdef _127MER_
struct kmer_t2 //use union later
{
	uint64_t kmer[4];
};
#endif



struct edge_node   // kmer-edge  the connection between sparse-kmer
{
	uint64_t edge: 50, edge_cov: 7, len: 6, used: 1;
	struct edge_node * nxt_edge;
};


//#pragma pack(4)    // do pack(4)  later
struct kmer_info
{
	//uint8_t used:1,split_left:1,split_right:1,removed:1,flip:1,marked:1,repeat:1;

	uint64_t used: 1, linear: 1, deleted: 1, single: 1, inEdge: 2, twin: 2, cov1: 16; //added for soapdenovo

	//uint16_t cov1:16;

	//uint32_t edge_id;//added for soapdenovo

	struct edge_node * left;
	struct edge_node * right;



};
struct kmer_info_r1
{
	uint16_t cov1: 16;
};

struct bucket2      //sparse-kmer
{
	struct kmer_t2 kmer_t2;
	struct kmer_info kmer_info;
	bucket2 * nxt_bucket;
};

struct bucket2_r1  //sparse-kmer struct  for round1 ,
{
	struct kmer_t2 kmer_t2;
	struct kmer_info_r1 kmer_info;
	bucket2_r1 * nxt_bucket;
};

struct hashtable2
{
	struct bucket2 ** store_pos;
	size_t ht_sz;
};

struct read_t //reads bits struct ...
{
	uint64_t read_bits[100];
	int readLen;
};



//function for hashtable
//void Init_HT2(struct hashtable2* ht,size_t ht_sz);
//bool look_up_in_a_list2_r1(struct kmer_t2 *seq,struct bucket2_r1 *** ptr);
//bool look_up_in_a_list2(struct kmer_t2 *seq,struct bucket2 *** ptr);
//void free_hashtable(hashtable2 *ht);

inline void Init_HT2 ( struct hashtable2 * ht, size_t ht_sz )
{
	ht->ht_sz = ht_sz;
	ht->store_pos = ( struct bucket2 ** ) calloc ( ht_sz, sizeof ( struct bucket2 * ) );

	for ( size_t i = 0; i < ht_sz; ++i )
	{
		ht->store_pos[i] = NULL;
	}
}


inline bool look_up_in_a_list2_r1 ( struct kmer_t2 * seq, struct bucket2_r1 ** * ptr )
{
	while ( ( **ptr ) != NULL )
	{
		if ( memcmp ( & ( ( **ptr )->kmer_t2.kmer ), & ( seq->kmer ), sizeof ( kmer_t2 ) ) == 0 )
		{
			break;
		}

		( *ptr ) = & ( ( **ptr )->nxt_bucket );
	}

	return ( ( **ptr ) != NULL );
}

inline bool look_up_in_a_list2 ( struct kmer_t2 * seq, struct bucket2 ** * ptr )
{
	while ( ( **ptr ) != NULL )
	{
		if ( memcmp ( & ( ( **ptr )->kmer_t2.kmer ), & ( seq->kmer ), sizeof ( kmer_t2 ) ) == 0 )
		{
			break;
		}

		( *ptr ) = & ( ( **ptr )->nxt_bucket );
	}

	return ( ( **ptr ) != NULL );
}

inline void free_bucket ( bucket2 * tmp )
{
	edge_node * edge, *edge2;
	edge = tmp->kmer_info.left;

	while ( edge )
	{
		edge2 = edge;
		edge = edge->nxt_edge;
		free ( edge );
	}

	edge = tmp->kmer_info.right;

	while ( edge )
	{
		edge2 = edge;
		edge = edge->nxt_edge;
		free ( edge );
	}

	free ( tmp );
}
inline void free_hashtable ( hashtable2 * ht )
{
	bucket2 * tmp, *tmp2;

	for ( size_t i = 0; i < ht->ht_sz; ++i )
	{
		tmp = ( ht->store_pos ) [i];

		while ( tmp )
		{
			tmp2 = tmp;
			tmp = tmp->nxt_bucket;
			free ( tmp2 );
		}
	}

	free ( ht->store_pos );
}




#endif

/*
 * build_edge.cpp
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

#include "core.h"
#include "stdinc.h"
#include "build_edge.h"
#include "seq_util.h"
#include "zlib.h"
#include "sparse_kmer.h"

/*************************************************
Function:
    RemovingWeakNodesAndEdges2
Description:
    1. Removes the nodes whose frequencies are not greater than threadhold NodeCovTh.
    2. Removes the kmer-edges (the connections between nodes) whose frequencies
    are not greater than threadhold EdgeCovTh.

    Some statistics function here temporarily added.
Input:
    1. ht:      the graph hash table
    2. K_size:      kmer size
    3. NodeCovTh:   the threadhold for removing low coverage nodes
    4. EdgeCovTh:   the threadhold for removing low coverage kmer-edges
    5. bucket_cnt:      the current bucket number
    6. edge_cnt:        the current kmer-edge number
Output:
    1. bucket_cnt:      the bucket number after removing weak nodes and kmer-edges
    2. edge_cnt:        the kmer-edge number after removing weak nodes and kmer-edges
Return:
    None.
*************************************************/
void RemovingWeakNodesAndEdges2 ( hashtable2 * ht, int K_size, int NodeCovTh, int EdgeCovTh, size_t * bucket_cnt, size_t * edge_cnt )
{
	stat_edge_num ( ht );
	stat_edge_cvg_len ( ht );
	int Removed_Nodes_cnt = 0, Removed_Edges_cnt = 0;
	bucket2 * bktptr = NULL, *bktptr_tmp = NULL;
	bucket2 ** bktp2p = NULL;
	edge_node * edge_ptr = NULL, *next_edge = NULL, *edge_tmp = NULL;
	int smaller;
	fprintf ( stderr, "Start to remove weak nodes and kmer-edges.\n" );

	/*
	for(size_t i=0;i<ht->ht_sz;++i)
	{
	    bktptr=ht->store_pos[i];
	    while(bktptr!=NULL)
	    {
	        if(bktptr->kmer_info.cov1==0)printf("zero\n");

	        bktptr=bktptr->nxt_bucket;
	    }

	}*/

	//removing weak nodes
	for ( size_t i = 0; i < ht->ht_sz; ++i )
	{
		bktptr = ht->store_pos[i];

		while ( bktptr != NULL )
		{
			if ( bktptr->kmer_info.cov1 <= NodeCovTh )
			{
				bktptr->kmer_info.deleted = 1;
				Removed_Nodes_cnt++;
				edge_ptr = bktptr->kmer_info.right;

				while ( edge_ptr )
				{
					edge_ptr->used = 1;
					edge_ptr = edge_ptr->nxt_edge;
					Removed_Edges_cnt++;
				}

				edge_ptr = bktptr->kmer_info.left;

				while ( edge_ptr )
				{
					edge_ptr->used = 1;
					edge_ptr = edge_ptr->nxt_edge;
					Removed_Edges_cnt++;
				}
			}

			bktptr = bktptr->nxt_bucket;
		}
	}

	//removing dead edges
	for ( size_t i = 0; i < ht->ht_sz; ++i )
	{
		bktptr = ht->store_pos[i];

		while ( bktptr != NULL )
		{
			edge_ptr = bktptr->kmer_info.right;

			while ( edge_ptr )
			{
				if ( edge_ptr->edge_cov <= EdgeCovTh )
				{
					edge_ptr->used = 1; //becasuse the cvg of edges is symmetrial, so it's ok
					Removed_Edges_cnt++;
				}
				else
				{
					bktptr_tmp = lastKmer ( ht, K_size, bktptr, edge_ptr, 0, smaller );

					if ( !bktptr_tmp )
					{
						fprintf ( stderr, "ERROR: to node not found error!\n" );
						exit ( -1 );
					}

					if ( bktptr_tmp ->kmer_info.deleted )
					{
						edge_ptr->used = 1;
						Removed_Edges_cnt++;
					}
				}

				edge_ptr = edge_ptr->nxt_edge;
			}

			edge_ptr = bktptr->kmer_info.left;

			while ( edge_ptr )
			{
				if ( edge_ptr->edge_cov <= EdgeCovTh )
				{
					edge_ptr->used = 1; //becasuse the cvg of edges is symmetrial, so it's ok
					Removed_Edges_cnt++;
				}
				else
				{
					bktptr_tmp = lastKmer ( ht, K_size, bktptr, edge_ptr, 1, smaller );

					if ( !bktptr_tmp )
					{
						fprintf ( stderr, "ERROR: to node not found error! \n" );
						exit ( -1 );
					}

					if ( bktptr_tmp ->kmer_info.deleted )
					{
						edge_ptr->used = 1;
						Removed_Edges_cnt++;
					}
				}

				edge_ptr = edge_ptr->nxt_edge;
			}

			bktptr = bktptr->nxt_bucket;
		}
	}

	for ( size_t i = 0; i < ht->ht_sz; ++i )
	{
		bktptr = ht->store_pos[i];
		bktp2p = & ( ht->store_pos[i] );

		while ( bktptr != NULL )
		{
			edge_ptr = bktptr->kmer_info.right;

			while ( edge_ptr )
			{
				next_edge = edge_ptr->nxt_edge;

				if ( edge_ptr->used )
				{
					removeEdge ( bktptr, edge_ptr, 0 );
					//Removed_Edges_cnt2++;
				}

				edge_ptr = next_edge;
			}

			edge_ptr = bktptr->kmer_info.left;

			while ( edge_ptr )
			{
				next_edge = edge_ptr->nxt_edge;

				if ( edge_ptr->used )
				{
					removeEdge ( bktptr, edge_ptr, 1 );
					//Removed_Edges_cnt2++;
				}

				edge_ptr = next_edge;
			}

			bktptr_tmp = bktptr->nxt_bucket;

			if ( bktptr->kmer_info.deleted )
			{
				free ( bktptr );
				( *bktp2p ) = bktptr_tmp;
				//Removed_Nodes_cnt2++;
			}
			else
			{
				bktp2p = & ( bktptr->nxt_bucket );
			}

			bktptr = bktptr_tmp;
		}
	}

	fprintf ( stderr, "%llu nodes removed.\n", Removed_Nodes_cnt );
	fprintf ( stderr, "%llu edges removed.\n", Removed_Edges_cnt );
	fprintf ( stderr, "\n" );
	( *bucket_cnt ) -= Removed_Nodes_cnt;
	( *edge_cnt ) -= Removed_Edges_cnt;
}



/*************************************************
Function:
    removeMinorTips
Description:
    Removes minor tips whose length are not longer than a threadhold usually set to 2*Kmer.
Input:
    1. ht:      the graph hash table
    2. K_size:      kmer size
    3. cut_len_tip:     the threadhold for tips' length
Output:
    1. tip_c:       useless
Return:
    None.
*************************************************/
void removeMinorTips ( struct hashtable2 * ht, int K_size, int cut_len_tip, int & tip_c )
{
	mask1in1out ( ht );
	bucket2 * bktptr = NULL;
	size_t flag = 1;
	size_t total = 0;
	int j = 0;

	while ( flag )
	{
		flag = 0;

		for ( size_t i = 0; i < ht->ht_sz; ++i )
		{
			bktptr = ht->store_pos[i];

			while ( bktptr != NULL )
			{
				flag += clipTipFromNode ( ht, K_size, bktptr, cut_len_tip );
				bktptr = bktptr->nxt_bucket;
			}
		}

		j++;

		if ( flag )
		{
			fprintf ( stderr, "%llu tips removed in cycle %d.\n\n", flag, j );
			total += flag;
		}
		else
		{
			fprintf ( stderr, "Total %llu tips removed.\n", total );
		}

		if ( flag ) { mask1in1out ( ht ); }
	}
}


/*************************************************
Function:
    mask1in1out
Description:
    Marks "1 in 1 out" node linear.
Input:
    1. ht:      the graph hashtable
Output:
    None.
Return:
    None.
*************************************************/
static void mask1in1out ( hashtable2 * ht )
{
	size_t total = 0, linear = 0;
	static int call_times;
	call_times++;

	for ( size_t i = 0; i < ht->ht_sz; ++i )
	{
		struct bucket2 * bkt_ptr = ht->store_pos[i];

		while ( bkt_ptr )
		{
			total++;//for stat

			if ( ( bkt_ptr->kmer_info.left != NULL && bkt_ptr->kmer_info.left->nxt_edge == NULL )
			        && ( bkt_ptr->kmer_info.right != NULL && bkt_ptr->kmer_info.right->nxt_edge == NULL ) )
			{
				bkt_ptr->kmer_info.linear = 1;
				linear++;//for stat
			}
			else
			{
				bkt_ptr->kmer_info.linear = 0;
			}

			bkt_ptr = bkt_ptr->nxt_bucket;
		}
	}

	//fprintf(stderr,"Masking linear nodes, times: %d\n",call_times);
	fprintf ( stderr, "Total nodes number: %llu\n", total );
	fprintf ( stderr, "Linear nodes number: %llu\n", linear );
}


/*************************************************
Function:
    clipTipFromNode
Description:
    Removes tips from an end node.
Input:
    1. ht:      the graph hashtable
    2. K_size:      the kmer size
    3. node:            the tips starting node
    4. cut_len_tip:     the threadhold of tips' length
Output:
    None.
Return:
    1 if clips a tip successfully.
*************************************************/
static int clipTipFromNode ( hashtable2 * ht, int K_size, bucket2 * node, int cut_len_tip ) //only for remove minor tips
{
	//linear return  0
	if ( node->kmer_info.linear || node->kmer_info.deleted )
	{
		return 0;
	}

	// for not linear
	int in_num, out_num;
	int sum_edge_len = 0;
	int smaller;
	bool is_left;
	bool pre_is_left;
	edge_node * edge0;
	in_num = count_left_edge_num ( node );
	out_num = count_right_edge_num ( node );

	if ( in_num == 0 && out_num == 1 ) { is_left = 0; }
	else if ( in_num == 1 && out_num == 0 ) { is_left = 1; }
	else { return 0; }

	if ( is_left ) { edge0 = node->kmer_info.left; }
	else { edge0 = node->kmer_info.right; }

	bucket2 * next, *pre_node;
	pre_node = node;
	next = lastKmer ( ht, K_size, node, edge0, is_left, smaller );

	while ( next->kmer_info.linear )
	{
		if ( sum_edge_len > cut_len_tip ) { return 0; }

		is_left = ! ( is_left ^ smaller );

		if ( is_left ) { edge0 = next->kmer_info.left; }
		else { edge0 = next->kmer_info.right; }

		sum_edge_len += edge0->len + 1;
		pre_node = next;
		next = lastKmer ( ht, K_size, next, edge0, is_left, smaller );

		if ( !next )
		{
			fprintf ( stderr, "ERROR: linear edge not found error !\n" );
			exit ( -1 );
		}
	}

	pre_is_left = is_left;
	is_left = ( is_left ^ smaller ); //back check orientation...
	in_num = count_left_edge_num ( next );
	out_num = count_right_edge_num ( next );

	if ( is_left ) //check the last node left branch or not
	{
		if ( in_num == 1 )
		{
			return 0;
		}
		else if ( in_num > 1 )
		{
			edge_node * edge1 = NULL, * temp_edge = NULL;
			bucket2 * temp_bucket = NULL;
			int max_cvg = 0, single_cvg = 0, temp_smaller;
			temp_edge = next->kmer_info.left;

			while ( temp_edge )
			{
				single_cvg = temp_edge->edge_cov;

				if ( single_cvg > max_cvg ) { max_cvg = single_cvg; }

				if ( !edge1 )
				{
					temp_bucket = lastKmer ( ht, K_size, next, temp_edge, 1, temp_smaller );

					if ( !temp_bucket )
					{
						fprintf ( stderr, "ERROR: edge to NULL found error ! a\n" );
						exit ( 1 );
					}

					if ( pre_node == temp_bucket ) { edge1 = temp_edge; }
				}

				temp_edge = temp_edge->nxt_edge;
			}

			if ( !edge1 )
			{
				fprintf ( stderr, "ERROR: edge to node not found error ! b\n" );
				exit ( 1 );
			}

			if ( edge1->edge_cov < max_cvg )
			{
				removeEdge ( next, edge1, 1 );
				removeEdge ( pre_node, edge0, pre_is_left );
				node->kmer_info.deleted = 1;
				pre_node->kmer_info.deleted = 1;
				return 1;
			}
			else
			{
				return 0;
			}
		}
		else
		{
			fprintf ( stderr, "ERROR: left tips oritation error or edge not found error ! a\n" );
			exit ( -1 );
		}
	}
	else
	{
		if ( out_num == 1 )
		{
			return 0;
		}
		else if ( out_num > 1 )
		{
			edge_node * edge1 = NULL, * temp_edge = NULL;
			bucket2 * temp_bucket = NULL;
			int max_cvg = 0, single_cvg = 0, temp_smaller;
			//ok change it to a edge_remove thred_hold locally later
			//or only if it is the least cvg ->remove it
			temp_edge = next->kmer_info.right;

			while ( temp_edge )
			{
				single_cvg = temp_edge->edge_cov;

				if ( single_cvg > max_cvg ) { max_cvg = single_cvg; }

				if ( !edge1 )
				{
					temp_bucket = lastKmer ( ht, K_size, next, temp_edge, 0, temp_smaller );

					if ( !temp_bucket )
					{
						fprintf ( stderr, "ERROR: edge to NULL found, error ! b\n" );
						exit ( -1 );
					}

					if ( pre_node == temp_bucket ) { edge1 = temp_edge; }
				}

				temp_edge = temp_edge->nxt_edge;
			}

			if ( !edge1 )
			{
				fprintf ( stderr, "ERROR: edge to node not found error ! e\n" );
				exit ( 1 );
			}

			if ( edge1->edge_cov < max_cvg )
			{
				removeEdge ( next, edge1, 0 );
				removeEdge ( pre_node, edge0, pre_is_left );
				node->kmer_info.deleted = 1;
				pre_node->kmer_info.deleted = 1;
				return 1;
			}
			else
			{
				return 0;
			}
		}
		else
		{
			fprintf ( stderr, "ERROR: right tips oritation error or edge not found error! b\n" );
			exit ( -1 );
		}
	}
}



static int count_left_edge_num ( bucket2 * bkt ) //63 127 same
{
	int ret = 0;

	if ( bkt )
	{
		edge_node * left_edge = bkt->kmer_info.left;

		while ( left_edge )
		{
			ret++;
			left_edge = left_edge->nxt_edge;
		}
	}

	return ret;
}

static int count_right_edge_num ( bucket2 * bkt ) //63 127 same
{
	int ret = 0;

	if ( bkt )
	{
		edge_node * right_edge = bkt->kmer_info.right;

		while ( right_edge )
		{
			ret++;
			right_edge = right_edge->nxt_edge;
		}
	}

	return ret;
}


/*************************************************
Function:
    lastKmer
Description:
    Searches the node that a node's kmer-edge end with.
Input:
    1. ht:      the graph hashtable
    2. K_size:      kmer size
    3. node:        the node whose kmer-edge will be searched
    4. edge:        the kmer-edge
    5. is_left:     whether the kmer-edge on the node's left side
Output:
    1. smaller:     whether the searched result, a kmer is smaller than its reversed complement
Return:
    A pointer to the found node.
    Null if not found.
*************************************************/
static bucket2 * lastKmer ( hashtable2 * ht, int K_size, bucket2 * node, edge_node * edge, int is_left, int & smaller ) //NEW
{
	if ( !node || !edge ) { return NULL; }

	kmer_t2 t_kmer, f_kmer;
	t_kmer = node->kmer_t2;
	kmer_t2 edge_seq;
	memset ( edge_seq.kmer, 0, sizeof ( edge_seq ) );
	( edge_seq.kmer ) [sizeof ( edge_seq ) / sizeof ( uint64_t ) - 1] = edge->edge;
	int edge_len = edge->len + 1;

	if ( edge_len > K_size )
	{
		fprintf ( stderr, "ERROR: g value should be no great than kmer size!\n" );
		exit ( -1 );
	}

	kmer_t2 KMER_FILTER;
	initKmerFilter ( K_size, &KMER_FILTER );

	if ( is_left ) //left edge
	{
		kmerMoveRight ( &t_kmer, edge_len );
		kmerMoveLeft ( &edge_seq, K_size - edge_len );
		kmerOr ( &t_kmer, &edge_seq );
		kmerAnd ( &t_kmer, &KMER_FILTER );
	}
	else
	{
		kmerMoveLeft ( &t_kmer, edge_len );
		kmerOr ( &t_kmer, &edge_seq );
		kmerAnd ( &t_kmer, &KMER_FILTER );
	}

	f_kmer = t_kmer;
	reverseCompKmer ( &f_kmer, K_size );

	if ( kmerCompare ( &t_kmer, &f_kmer ) > 0 )
	{
		t_kmer = f_kmer;
		smaller = 0;
	}
	else { smaller = 1; }

	return search_kmer ( ht, &t_kmer );
}

static bucket2 * search_kmer ( hashtable2 * ht, kmer_t2 * t_kmer )
{
	uint64_t hv = MurmurHash64A ( t_kmer, sizeof ( kmer_t2 ), 0 );
	size_t hash_idx = ( size_t ) ( hv % ht->ht_sz );
	bucket2 * starter = ht->store_pos[hash_idx];

	while ( starter )
	{
		if ( kmerCompare ( & ( starter->kmer_t2 ), t_kmer ) == 0 )
		{
			return starter;
		}

		starter = starter->nxt_bucket;
	}

	return NULL;
}




static void removeEdge ( bucket2 * node, edge_node * edge, int is_left ) // remove only one side ...   //63 127 same ...
{
	edge_node * pre_edge = NULL, *cur_edge = NULL, *nxt_edge = NULL;

	if ( !node || !edge )
	{
		return ;
	}

	if ( is_left )
	{
		cur_edge = node->kmer_info.left;

		if ( cur_edge == NULL )
		{
			return ;
		}

		if ( cur_edge == edge )
		{
			nxt_edge = cur_edge->nxt_edge;
			free ( cur_edge );
			cur_edge = NULL;
			node->kmer_info.left = nxt_edge;
			return ;
		}
	}
	else
	{
		cur_edge = node->kmer_info.right;

		if ( cur_edge == NULL )
		{
			return ;
		}

		if ( cur_edge == edge )
		{
			nxt_edge = cur_edge->nxt_edge;
			free ( cur_edge );
			cur_edge = NULL;
			node->kmer_info.right = nxt_edge;
			return ;
		}
	}

	pre_edge = cur_edge;
	cur_edge = cur_edge->nxt_edge;

	while ( cur_edge )
	{
		if ( cur_edge == edge ) { break; }

		pre_edge = cur_edge;
		cur_edge = cur_edge->nxt_edge;
	}

	if ( cur_edge )
	{
		nxt_edge = cur_edge->nxt_edge;
		free ( cur_edge );
		cur_edge = NULL;
		pre_edge->nxt_edge = nxt_edge;
	}
}



static void stat_edge_num ( hashtable2 * ht ) //63 127 same
{
	int l_num = 0, r_num = 0;
	size_t total_edge_num = 0, total_node_num = 0;
	bucket2 * bkt = NULL;
	map<int, size_t> edge_num_map;

	for ( size_t i = 0; i < ht->ht_sz; i++ )
	{
		bkt = ht->store_pos[i];

		while ( bkt )
		{
			total_node_num++;
			l_num = count_left_edge_num ( bkt );
			r_num = count_right_edge_num ( bkt );
			total_edge_num += ( l_num + r_num );
			edge_num_map[l_num]++;
			edge_num_map[r_num]++;
			bkt = bkt->nxt_bucket;
		}
	}

	ofstream o_edge_num ( "edge_num_stat.txt" );
	o_edge_num << "Total nodes number:" << total_node_num << endl;
	o_edge_num << "Total kmer-edges number:" << total_edge_num << endl;
	o_edge_num << "Average kmer-edges number per node:" << ( double ) total_edge_num / total_node_num << endl;
	o_edge_num << "The frequence of kmer-edges number on a node's one side as below :" << endl;
	map<int, size_t>::iterator it;

	for ( it = edge_num_map.begin(); it != edge_num_map.end(); ++it )
	{
		o_edge_num << it->first << "\t" << it->second << endl;
	}

	o_edge_num.close();
}



static void stat_edge_cvg_len ( hashtable2 * ht )
{
	map<int, size_t> edge_cvg_map;
	map<int, size_t> edge_len_map;
	bucket2 * bkt = NULL;
	edge_node * temp_edge = NULL;

	for ( size_t i = 0; i < ht->ht_sz; i++ )
	{
		bkt = ht->store_pos[i];

		while ( bkt )
		{
			//left
			temp_edge = bkt->kmer_info.left;

			while ( temp_edge )
			{
				edge_cvg_map[temp_edge->edge_cov]++;
				edge_len_map[temp_edge->len]++;
				temp_edge = temp_edge->nxt_edge;
			}

			//right
			temp_edge = bkt->kmer_info.right;

			while ( temp_edge )
			{
				edge_cvg_map[temp_edge->edge_cov]++;
				edge_len_map[temp_edge->len]++;
				temp_edge = temp_edge->nxt_edge;
			}

			bkt = bkt->nxt_bucket;
		}
	}

	ofstream o_edge_cvg ( "edge_cvg_stat.txt" );
	ofstream o_edge_len ( "edge_len_stat.txt" );
	map<int, size_t>::iterator it;

	for ( it = edge_cvg_map.begin(); it != edge_cvg_map.end(); ++it )
	{
		o_edge_cvg << it->first << "\t" << it->second << endl;
	}

	for ( it = edge_len_map.begin(); it != edge_len_map.end(); ++it )
	{
		o_edge_len << it->first << "\t" << it->second << endl;
	}

	o_edge_cvg.close();
	o_edge_len.close();
}



/*************************************************
Function:
    kmer2edges
Description:
    This is the main function for building edges by compacting the linear nodes.
Input:
    1. ht:      the graph hashtable
    2. K_size:      kmer size
    3. out_file:    the name of output file containing edges sequence
Output:
    None.
Return:
    None.
*************************************************/
void kmer2edges ( hashtable2 * ht, int K_size, char * outfile )
{
	FILE * fp;
	char temp[256];
	sprintf ( temp, "%s", outfile );
	fp = fopen ( temp, "w" );

	if ( fp == NULL )
	{
		fprintf ( stderr, "ERROR: Can't create file %s. \n", temp );
		exit ( -1 );
	}

	make_edge ( ht, K_size, fp );
	fclose ( fp );
}

static void make_edge ( hashtable2 * ht, int K_size, FILE * fp ) //63 127 same
{
	bucket2 * bktptr;

	for ( size_t i = 0; i < ht->ht_sz; ++i )
	{
		bktptr = ht->store_pos[i];

		while ( bktptr != NULL )
		{
			startEdgeFromNode ( ht, K_size, bktptr, fp );
			bktptr = bktptr->nxt_bucket;
		}
	}
}


/*************************************************
Function:
    startEdgeFromNode
Description:
    Constructs edges from a branched node or end node.
    for every branch (left , right)
    1. Puts the linear node into a stack
    2. Checks the edge to be built form the stack are plalindrome or not
    3. Builds an edge by merge the linear nodes
Input:
    1. ht:      the graph hashtable
    2. K_size:      kmer size
    3. fp:      the file pointer for writing out edge sequences
Output:
    None.
Return:
    Zero.
*************************************************/
static int startEdgeFromNode ( hashtable2 * ht, int K_size, bucket2 * node, FILE * fp )
{
	static size_t call_times;
	call_times++;

	if ( node->kmer_info.linear || node->kmer_info.deleted )
	{
		return 0;//linear node ...
	}

	int left, right;
	left = count_left_edge_num ( node );
	right = count_right_edge_num ( node );

	if ( left == 0 && right == 0 )
	{
		return 0; //it's a dead node
	}

	list<stacked_node2 *> stack;
	edge_node * t_edge = NULL, *t_next = NULL;
	stacked_node2 * t_stacked_node = NULL;
	vector<preEDGE2> loops_edges;
	int node_c;
	//for right edge
	t_edge = node->kmer_info.right;

	while ( t_edge )
	{
		if ( t_edge->used == 1 )
		{
			t_edge = t_edge->nxt_edge;
			continue;
		}

		t_stacked_node = ( stacked_node2 * ) malloc ( sizeof ( stacked_node2 ) );
		t_stacked_node->node = node;
		t_stacked_node->is_left = 0;
		t_stacked_node->edge = t_edge;
		t_stacked_node->next = NULL;
		stack.push_back ( t_stacked_node );
		t_edge->used = 1;
		stringBeads ( ht, K_size, stack, t_stacked_node, t_edge, &node_c );
		process_1stack ( ht, K_size, stack, fp, loops_edges );
		t_next = t_edge->nxt_edge;//because this procedure will remove the edge t_edge
		dislink ( ht, K_size, stack.front() );

		if ( stack.size() > 2 )
		{
			stack.pop_back();//change the stack

			if ( stack.back() && stack.size() > 1 ) //last but second node
			{
				dislink ( ht, K_size, stack.back() );
			}
		}

		stacked_node2 * head, *tmp_node;
		head = stack.front();

		while ( head )
		{
			tmp_node = head;
			free ( tmp_node );
			head = head->next;
		}

		stack.clear();
		t_edge = t_next;
	}

	//for left edge
	t_edge = node->kmer_info.left;

	while ( t_edge )
	{
		if ( t_edge->used == 1 )
		{
			t_edge = t_edge->nxt_edge;
			continue;
		}

		t_stacked_node = ( stacked_node2 * ) malloc ( sizeof ( stacked_node2 ) );
		t_stacked_node->node = node;
		t_stacked_node->is_left = 1;
		t_stacked_node->edge = t_edge;
		t_stacked_node->next = NULL;
		stack.push_back ( t_stacked_node );
		t_edge->used = 1;
		stringBeads ( ht, K_size, stack, t_stacked_node, t_edge, &node_c ); //
		process_1stack ( ht, K_size, stack, fp, loops_edges );
		t_next = t_edge->nxt_edge;//because this procedure will remove the edge t_edge
		dislink ( ht, K_size, stack.front() );

		if ( stack.size() > 2 )
		{
			stack.pop_back();//change the stack

			if ( stack.back() && stack.size() > 1 ) //last but second node
			{
				dislink ( ht, K_size, stack.back() );
			}
		}

		//debug<<"before free stack"<<endl;
		stacked_node2 * head, *tmp_node;
		head = stack.front();

		while ( head )
		{
			tmp_node = head;
			free ( tmp_node );
			head = head->next;
		}

		stack.clear();
		t_edge = t_next;
	}

	if ( loops_edges.size() > 0 )
	{
		//fprintf(stderr,"loops_edges size %llu\n",loops_edges.size());
		int i, j, size;
		bool need_output;
		size = loops_edges.size();
		need_output = 1;

		//bool debug = 0;
		for ( i = 0; i < size; i++ )
		{
			string seq = * ( loops_edges[i].full_edge );
			string rc_seq = revCompSeq ( seq );
			/*
			if(seq.compare("AATTGGACGTGAGAGCAAATTGTATTGAGCATACAATTTGCTCTCACGTCCAATT") == 0) {
			                    fprintf(stderr,"in loops_edges %d %s\n",i,seq.c_str());
			    debug = 1;
			            }

			            if(seq.compare("AATTGGACGTGAGAGCAAATTGTATGCTCAATACAATTTGCTCTCACGTCCAATT") == 0) {
			                    fprintf(stderr,"in loops_edges %d %s\n",i,seq.c_str());
			    debug = 1;
			            }

			if(debug ){
			    fprintf(stderr, "%d %s\n",i,seq.c_str());
			    fprintf(stderr, "%d %s\n",i,rc_seq.c_str());
			}*/

			for ( j = i + 1; j < size; j++ )
			{
				string cur_seq = * ( loops_edges[j].full_edge );

				if ( seq.compare ( cur_seq ) == 0 )
				{
					fprintf ( stderr, "ERROR: two equal loop edge sequence from same node, this should not happen!\n" );
					fprintf ( stderr, "%s\n", seq.c_str() );
					exit ( -1 );
				}

				if ( rc_seq.compare ( cur_seq ) == 0 )
				{
					fprintf ( stderr, "INFO: two loop edge sequence are reversed complemental!\n" );
					fprintf ( stderr, "%s\n", seq.c_str() );
					fprintf ( stderr, "%s\n", rc_seq.c_str() );
					need_output = 0;
					loops_edges[j].cvg += loops_edges[i].cvg;
					break;
				}
			}

			if ( need_output )
			{
				output_1edge ( &loops_edges[i], K_size, fp );
				//fprintf(stderr,"need output %d %s\n",i,seq.c_str());
			}

			delete ( loops_edges[i].full_edge );
			need_output = 1;
		}
	}

	return 0;
}



/*************************************************
Function:
    stringBeads
Description:
    Puts the linear nodes into a stack.
Input:
    1. ht:      the graph hashtalbe
    2. K_size:      kmer size
    3. stack:       a stack
    4. from_node:       the starting node of the stack
    5. from_edge:       the kmer-edge of the first node
    6. node_c:      useless
Output:
    None.
Return:
    None.
*************************************************/
static void stringBeads ( hashtable2 * ht, int K_size, list<stacked_node2 *> &stack, stacked_node2 * from_node, edge_node * from_edge, int * node_c )
{
	static size_t call_times;
	call_times++;
	bucket2 * t_bucket = from_node->node;
	edge_node * t_edge = from_edge;
	stacked_node2 * t_stacked_node = from_node;
	int is_left = from_node->is_left;
	int t_smaller;
	t_edge->used = 1;
	t_bucket = lastKmer ( ht, K_size, t_bucket, t_edge, is_left, t_smaller );

	if ( !t_bucket )
	{
		fprintf ( stderr, "ERROR: to node not found in stringBeads()\n" );
		exit ( -1 );
	}

	while ( t_bucket && t_bucket->kmer_info.linear )
	{
		t_stacked_node = ( stacked_node2 * ) malloc ( sizeof ( stacked_node2 ) );
		t_stacked_node->node = t_bucket;
		is_left = ! ( is_left ^ t_smaller );
		t_stacked_node->is_left = is_left;

		if ( is_left ) { t_stacked_node->edge = t_bucket->kmer_info.left; }
		else { t_stacked_node->edge = t_bucket->kmer_info.right; }

		t_stacked_node->next = NULL;
		( ( stacked_node2 * ) stack.back() )->next = t_stacked_node;
		stack.push_back ( t_stacked_node );
		t_stacked_node->edge->used = 1;
		t_bucket = lastKmer ( ht, K_size, t_bucket, t_stacked_node->edge, is_left, t_smaller );
	}

	if ( t_bucket ) //should be always true   for end node ..
	{
		t_stacked_node = ( stacked_node2 * ) malloc ( sizeof ( stacked_node2 ) );
		t_stacked_node->node = t_bucket;
		is_left = ! ( is_left ^ t_smaller );
		t_stacked_node->is_left = is_left;
		t_stacked_node->edge = NULL;
		t_stacked_node->next = NULL;
		( ( stacked_node2 * ) stack.back() )->next = t_stacked_node;
		stack.push_back ( t_stacked_node );
	}
}
//for debug
static void pirntStack ( list<stacked_node2 *> &stack )
{
	static int times = 0;
	fprintf ( stderr, "call times %d \n ", times++ );
	stacked_node2 * ptr = stack.front();

	while ( ptr )
	{
		printKmer ( & ( ptr->node->kmer_t2 ), stderr );

		if ( ptr->edge )
			{ fprintf ( stderr, "%llx , %d ,", ptr->edge->edge, ptr->is_left ); }

		fprintf ( stderr, "->" );
		ptr = ptr->next;
	}

	fprintf ( stderr, "\n" );
}
/*************************************************
Function:
    process_1stack
Description:
    Processes the nodes in one stack
    1. Compacts the nodes to an edge
    2. Checks palindrome
    3. Calculates coverage
Input:
    1. ht:      the graph hashtable
    2. K_size:      kmer size
    3. stack:       the stack
    4. fp:      the file pointer for writing
Output:
    None.
Return:
    None.
*************************************************/
static void process_1stack ( hashtable2 * ht, int K_size, list<stacked_node2 *> &stack, FILE * fp, vector<preEDGE2> &loops_edges )
{
	static size_t edge_c;// edge id
	static preEDGE2 long_edge_buf;
	preEDGE2 loops;
	int TipLenTh = 3 * K_size; //orig 100
	int TipCovTh = 5;

	if ( stack.size() < 2 )
	{
		fprintf ( stderr, "only %llu nodes in the stack \n", stack.size() );
		exit ( -1 );
	}
	else
	{
		//palindrome check
		string full_edge = stack2string ( ht, K_size, stack ); //when output  skip the first kmer first
		stacked_node2 * test = stack.front();
		bool palindrome = check_palindrome ( full_edge );
		int bal_edge = !palindrome;
		stacked_node2 * from_node = stack.front();
		stacked_node2 * to_node = stack.back();
		long_edge_buf.from_node = from_node;
		long_edge_buf.to_node = to_node;
		long_edge_buf.full_edge = &full_edge;
		long_edge_buf.bal_edge = bal_edge;
		uint64_t symbol = 0; //cvg stat
		edge_c++;

		if ( stack.size() == 2 )
		{
			long_edge_buf.cvg = from_node->edge->edge_cov;
		}
		else
		{
			stacked_node2 * nd_tmp = from_node;

			while ( nd_tmp && nd_tmp->edge )
			{
				symbol += nd_tmp->edge->edge_cov * ( nd_tmp->edge->len + 1 );
				nd_tmp = nd_tmp->next;
			}

			int cvg = symbol / ( full_edge.size() - K_size );
			long_edge_buf.cvg = cvg;
		}

		int from_left, from_right, to_left, to_right;
		from_left = count_left_edge_num ( from_node->node );
		from_right = count_right_edge_num ( from_node->node );
		to_left = count_left_edge_num ( to_node->node );
		to_right = count_right_edge_num ( to_node->node );

		//tips control

		if ( ( ( from_left + from_right == 1 ) && ( to_left + to_right == 1 ) && ( full_edge.size() < TipLenTh ) )
		        || ( ( ( from_left + from_right == 1 ) || ( to_left + to_right == 1 ) )
		             && ( full_edge.size() < TipLenTh ) && long_edge_buf.cvg < TipCovTh ) ) //tips args
		{
			//if(full_edge.size()<TipLenTh && long_edge_buf.cvg<TipCovTh){//it's a tip or low cvg link
			static size_t tip_num;
			tip_num++;
		}
		else
		{
			//debug begin
			/*
			string bug_seq = *(long_edge_buf.full_edge);
			if(bug_seq.compare("AATTGGACGTGAGAGCAAATTGTATTGAGCATACAATTTGCTCTCACGTCCAATT") == 0) {
			    fprintf(stderr,"%s\n",bug_seq.c_str());
			    fprintf(stderr,"from %llx to %llx \n",long_edge_buf.from_node->node,long_edge_buf.to_node->node);

			}

			if(bug_seq.compare("AATTGGACGTGAGAGCAAATTGTATGCTCAATACAATTTGCTCTCACGTCCAATT") == 0) {
			                    fprintf(stderr,"%s\n",bug_seq.c_str());
			                    fprintf(stderr,"from %llx to %llx \n",long_edge_buf.from_node->node,long_edge_buf.to_node->node);

			            }*/

			//debug end
			if ( long_edge_buf.from_node->node == long_edge_buf.to_node->node )
			{
				loops = long_edge_buf;
				loops.full_edge = new string ( * ( long_edge_buf.full_edge ) );
				loops_edges.push_back ( loops );
			}
			else
			{
				//output edge
				output_1edge ( &long_edge_buf, K_size, fp );
			}
		}

		edge_c += bal_edge;
	}
}



// WARNING: the kmer atcg is different from soapdenovo's represent
static void output_1edge ( preEDGE2 * long_edge, int K_size, FILE * fp )
{
	fprintf ( fp, ">length %d,", long_edge->full_edge->size() - K_size );
	const char * seq = long_edge->full_edge->c_str();
	//uint64_t from_kmer[2],to_kmer[2];
	kmer_t2 from_kmer, to_kmer;
	get_kmer_from_seq ( seq, long_edge->full_edge->size() , K_size, 0, &from_kmer );
	get_kmer_from_seq ( seq, long_edge->full_edge->size() , K_size, long_edge->full_edge->size() - K_size, &to_kmer );
	uint64_t * from, *to;
	from = from_kmer.kmer;
	to = to_kmer.kmer;
#ifdef _63MER_
	fprintf ( fp, "%llx %llx,", from[0], from[1] );
	fprintf ( fp, "%llx %llx,", to[0], to[1] );
#endif
#ifdef _127MER_
	fprintf ( fp, "%llx %llx %llx %llx,", from[0], from[1], from[2], from[3] );
	fprintf ( fp, "%llx %llx %llx %llx,", to[0], to[1], to[2], to[3] );
#endif
	fprintf ( fp, "cvg %d,%d\n", long_edge->cvg, long_edge->bal_edge );
	fprintf ( fp, "%s", seq );
	fprintf ( fp, "\n" );
}


/*************************************************
Function:
    dislink
Description:
    Marks the kmer-edges between "form_node" and "form_node->next" used.
Input:
    1. ht:      the graph hashtable
    2. K_size:      kmer size
    3. from_node:       the first node of a stack
Output:
    None.
Return:
    None.
*************************************************/
static void dislink ( hashtable2 * ht, int K_size, stacked_node2 * from_node )
{
	from_node->edge->used = 1;
	int from_edge_len = from_node->edge->len;
	int smaller;

	if ( from_node->next )
	{
		stacked_node2 * from_next = from_node->next;

		if ( from_next->node == from_node->node && from_next->is_left != from_node->is_left )
		{
			return ;
		}

		if ( from_next->is_left ) //remove right edge
		{
			if ( from_next->node->kmer_info.linear )
			{
				bucket2 * node_tmp = lastKmer ( ht, K_size, from_next->node,  from_next->node->kmer_info.right,  0, smaller );

				if ( node_tmp == from_node->node )
				{
					from_next->node->kmer_info.right->used = 1;
				}
				else
				{
					fprintf ( stderr, "ERROR: to node not found in dislink()\n" );
				}
			}
			else
			{
				edge_node * edge_tmp = from_next->node->kmer_info.right;

				while ( edge_tmp )
				{
					bucket2 * node_tmp = lastKmer ( ht, K_size, from_next->node, edge_tmp,  0, smaller );

					if ( node_tmp == from_node->node && edge_tmp->len == from_edge_len ) //there may be two or more edges between two nodes
					{
						edge_tmp->used = 1;
						break;
					}

					edge_tmp = edge_tmp->nxt_edge;
				}

				if ( !edge_tmp )
				{
					fprintf ( stderr, "ERROR: to node not found in dislink()\n" );
				}
			}
		}
		else  // remove left edge
		{
			if ( from_next->node->kmer_info.linear )
			{
				bucket2 * node_tmp = lastKmer ( ht, K_size, from_next->node,  from_next->node->kmer_info.left,  1, smaller );

				if ( node_tmp == from_node->node )
				{
					from_next->node->kmer_info.left ->used = 1;
				}
				else
				{
					fprintf ( stderr, "ERROR: to node not found in dislink()\n" );
				}
			}
			else
			{
				edge_node * edge_tmp = from_next->node->kmer_info.left;

				while ( edge_tmp )
				{
					bucket2 * node_tmp = lastKmer ( ht, K_size, from_next->node, edge_tmp, 1, smaller );

					if ( node_tmp == from_node->node && edge_tmp->len == from_edge_len )
					{
						edge_tmp->used = 1;
						break;
					}

					edge_tmp = edge_tmp->nxt_edge;
				}

				if ( !edge_tmp )
				{
					fprintf ( stderr, "ERROR: to node not found in dislink()\n" );
				}
			}
		}
	}
}


/*************************************************
Function:
    stack2string
Description:
    Compacts the nodes in the stack to a sequence.
Input:
    1. ht:      the graph hashtable
    2. K_size:      kmer size
    3. stack:       the stack for processing
Output:
    None.
Return:
    The compacted string, namely the edge sequence.
*************************************************/          //63 127 differ, fixed
static string stack2string ( hashtable2 * ht, int K_size, list<stacked_node2 *> & stack )
{
	static size_t call_times;
	call_times++;
	string full_edge;
	stacked_node2 * t_stack_node = stack.front();
	char tmp[1024];
	uint64_t bits[2];
	kmer_t2 tmp_kmer = ( t_stack_node->node->kmer_t2 );

	if ( t_stack_node->is_left )
	{
		reverseCompKmer ( &tmp_kmer, K_size );
	}
	else
	{
	}

	bitsarr2str ( tmp_kmer.kmer, K_size, tmp, sizeof ( kmer_t2 ) / sizeof ( uint64_t ) );
	full_edge.append ( tmp ); //put first node

	while ( t_stack_node )
	{
		if ( t_stack_node->edge )
		{
			if ( t_stack_node->is_left )
			{
				bits[0] = get_rev_comp_seq ( t_stack_node->edge->edge, t_stack_node->edge->len + 1 );
				bitsarr2str ( bits, t_stack_node->edge->len + 1, tmp, 1 );
				full_edge.append ( tmp );
			}
			else
			{
				bits[0] = t_stack_node->edge->edge;
				bitsarr2str ( bits, t_stack_node->edge->len + 1, tmp, 1 );
				full_edge.append ( tmp );
			}
		}

		t_stack_node = t_stack_node->next;
	}

	return full_edge;
}

static bool check_palindrome ( string & str ) //63 127 same
{
	size_t size = str.size();
	size_t mid = ( size / 2 ) + 1;

	for ( size_t i = 0; i < mid; ++i )
	{
		switch ( str[i] )
		{
			case 'A':

				if ( str[size - i - 1] != 'T' ) { return 0; }

				break;
			case 'C':

				if ( str[size - i - 1] != 'G' ) { return 0; }

				break;
			case 'T':

				if ( str[size - i - 1] != 'A' ) { return 0; }

				break;
			case 'G':

				if ( str[size - i - 1] != 'C' ) { return 0; }

				break;
		}
	}

	return 1;
}

static string revCompSeq ( const string & str )
{
	string rc_seq;
	size_t size = str.size();

	for ( int i = size - 1; i >= 0; i-- )
	{
		switch ( str[i] )
		{
			case 'A':
				rc_seq.push_back ( 'T' );
				break;
			case 'C':
				rc_seq.push_back ( 'G' );
				break;
			case 'T':
				rc_seq.push_back ( 'A' );
				break;
			case 'G':
				rc_seq.push_back ( 'C' );
				break;
		}
	}

	return rc_seq;
}




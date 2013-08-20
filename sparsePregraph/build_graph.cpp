/*
 * build_graph.cpp
 *
 * This file is part of SOAPdenovo.
 *
 * Part of this file is refered and modified from SparseAssembler
 * (See <http://sourceforge.net/projects/sparseassembler/>).
 *
 */

#include "stdinc.h"
#include "core.h"
#include "global.h"

#include "seq_util.h"
#include "io_func.h"
#include "build_graph.h"

static void process_round1_threaded ( struct read_t * read, struct hashtable2 * ht, pthread_spinlock_t * locks, size_t * bucket_count, int K_size, int gap );
//static void process_round2_threaded(struct read_t *read,struct hashtable2 *ht,pthread_spinlock_t *locks,size_t *edge_cnt,int K_size,int gap);

//for debug
//static void process_round1_threaded_d(struct read_t *read, struct hashtable2 *ht,pthread_spinlock_t *locks,size_t *bucket_count,int K_size,int gap);
static void process_round2_threaded_d ( struct read_t * read, struct hashtable2 * ht, pthread_spinlock_t * locks, size_t * edge_cnt, int K_size, int gap );


/*************************************************
Function:
    run_process_threaded
Description:
    Builds the sparse de-Brujin graph from reads.
    It calls function "process_round1_threaded" for round1 and
    "process_round2_threaded_d" for round 2 building process.
Input:
    1. ht:      the graph hashtable
    2. locks:       the locks array for hashtable
    3. K_size:  kmer size
    3. gap:     the skipped distance
    4. read_num:        the number of reads for processing
    5. thrd_num_s:        the thread number for building sparse de-Brujin graph
    6. thrd_id:     current thread's id
    7. round:           current building round (1,2)
Output:
    None.
Return:
    None.
*************************************************/
void run_process_threaded ( struct hashtable2 * ht, pthread_spinlock_t * locks, int K_size, int gap, size_t read_num, int thrd_num_s, int thrd_id, int round )
{
	read_t read_tmp;

	for ( int i = thrd_id; i < read_num; i += thrd_num_s )
	{
		int bad_flag = 0;
		filter_N ( seq_t[i], bad_flag );

		if ( bad_flag ) {seq_t[i].clear(); continue;}

		Init_Read ( seq_t[i], read_tmp );

		if ( round == 1 )
		{
			//cout << "round 1:"<< seq_t[i] <<endl;
			process_round1_threaded ( &read_tmp, ht, locks, &bucket_count_total[thrd_id], K_size, gap );
		}
		else if ( round == 2 )
		{
			//cout << "round 2:"<< seq_t[i] <<endl;
			//process_round1_threaded_d(&read_tmp,ht,locks,&bucket_count_total[thrd_id],K_size,gap);
			process_round2_threaded_d ( &read_tmp, ht, locks, &edge_cnt_total[thrd_id], K_size, gap );
		}
		else
		{
			fprintf ( stderr, "ERROR: invalid round number!\n" );
			exit ( -1 );
		}
	}
}


/*************************************************
Function:
    process_round1_threaded
Description:
    Processes one read in round1:
    1. Chops a read to kmers
    2. Searches hash and select sparse-kmers
Input:
    1. read:        a read
    2. ht:      the graph hashtable
    3. locks:       the locks array for hashtable
    4. bucket_count:    useless
    5. K_size:      kmer size
    6. gap:         the skipped distance
Output:
    None.
Return:
    None.
*************************************************/
static void process_round1_threaded ( struct read_t * read, struct hashtable2 * ht, pthread_spinlock_t * locks, size_t * bucket_count, int K_size, int gap )
{
	int readLen = read->readLen;
	int OverlappingKmers = readLen - K_size + 1;

	if ( gap >= OverlappingKmers )
		{   return;}

	int Read_arr_sz = readLen / 32 + 1;
	int rem = readLen % 32;

	if ( rem == 0 )
		{Read_arr_sz--;}

#ifdef _63MER_
	int Kmer_arr_sz = 2;
	int tot_bits = Read_arr_sz * 64;
#endif
#ifdef _127MER_
	int Kmer_arr_sz = 4;
	int tot_bits = Read_arr_sz * 128;
#endif
	size_t ht_sz = ht->ht_sz;
	bool flip[500], found[500];
	size_t hash_idx[500];
	memset ( flip, 0, sizeof ( flip ) );
	memset ( found, 0, sizeof ( found ) );
	kmer_t2 seq[500], f_seq[500];
	memset ( seq, 0, sizeof ( seq ) );
	uint64_t hv[500], temp_bits[500];
	bucket2 ** bktptr[500];
	char c_str[500];

	for ( int j = 0; j < OverlappingKmers; j++ )
	{
		get_sub_arr ( read->read_bits, read->readLen, j, K_size, seq[j].kmer );
#ifdef _63MER_

		if ( K_size <= 31 ) //fix the represent bug
		{
			( seq[j].kmer ) [1] = ( seq[j].kmer ) [0];
			( seq[j].kmer ) [0] = 0;
		}

#endif
#ifdef _127MER_ //fix the represent bug

		if ( K_size <= 31 ) //fix the represent bug
		{
			( seq[j].kmer ) [3] = ( seq[j].kmer ) [0];
			( seq[j].kmer ) [0] = 0;
		}
		else if ( K_size <= 63 )
		{
			( seq[j].kmer ) [3] = ( seq[j].kmer ) [1];
			( seq[j].kmer ) [2] = ( seq[j].kmer ) [0];
			( seq[j].kmer ) [1] = 0;
			( seq[j].kmer ) [0] = 0;
		}
		else if ( K_size <= 95 )
		{
			( seq[j].kmer ) [3] = ( seq[j].kmer ) [2];
			( seq[j].kmer ) [2] = ( seq[j].kmer ) [1];
			( seq[j].kmer ) [1] = ( seq[j].kmer ) [0];
			( seq[j].kmer ) [0] = 0;
		}

#endif
		memcpy ( &f_seq[j], &seq[j], Kmer_arr_sz * sizeof ( uint64_t ) );
		get_rev_comp_seq_arr ( ( f_seq[j].kmer ), K_size, Kmer_arr_sz ); //TODO ,add 127mer support

		if ( uint64_t_cmp ( seq[j].kmer, f_seq[j].kmer, Kmer_arr_sz ) > 0 )
		{
			flip[j] = 1;
		}

		if ( flip[j] == 1 )
		{
			memcpy ( temp_bits, & ( seq[j].kmer ), Kmer_arr_sz * sizeof ( uint64_t ) );
			memcpy ( & ( seq[j].kmer ), & ( f_seq[j].kmer ), Kmer_arr_sz * sizeof ( uint64_t ) );
			memcpy ( & ( f_seq[j].kmer ), temp_bits, Kmer_arr_sz * sizeof ( uint64_t ) );
		}

		hv[j] = MurmurHash64A ( ( seq[j].kmer ), sizeof ( seq[j] ), 0 );
		hash_idx[j] = ( size_t ) ( hv[j] % ht_sz );
		bktptr[j] = & ( ht->store_pos[hash_idx[j]] );
	}

	int g, h;
	g = 0;

	for ( int k = 0; k < gap; ++k )
	{
		pthread_spin_lock ( &locks[hash_idx[k]] );
		found[k] = look_up_in_a_list2_r1 ( &seq[k], ( struct bucket2_r1 ** * ) &bktptr[k] );
		pthread_spin_unlock ( &locks[hash_idx[k]] );

		if ( found[k] == 1 )
		{
			g = k;
			break;
		}
	}

	for ( int j = g; j < OverlappingKmers; )
	{
		h = gap;

		for ( int k = 0; k < gap; ++k )
		{
			if ( ( j + k ) >= OverlappingKmers - 1 )
			{
				h = k + 1; //舍弃掉最后一个kmer
				break;
			}

			pthread_spin_lock ( &locks[hash_idx[j + k]] );
			found[j + k] = look_up_in_a_list2_r1 ( &seq[j + k], ( bucket2_r1 ** * ) &bktptr[j + k] ); //lock...
			pthread_spin_unlock ( &locks[hash_idx[j + k]] );

			if ( k > 0 && found[j + k] == 1 )
			{
				h = k;
				break;
			}
		}

		pthread_spin_lock ( &locks[hash_idx[j]] );
		found[j] = look_up_in_a_list2_r1 ( &seq[j], ( bucket2_r1 ** * ) &bktptr[j] ); //lock...
		pthread_spin_unlock ( &locks[hash_idx[j]] );

		if ( found[j] == 0 )
		{
			pthread_spin_lock ( &locks[hash_idx[j]] );
			* ( bktptr[j] ) = ( struct bucket2 * ) malloc ( sizeof ( struct bucket2_r1 ) ); //lock ...
			memset ( * ( bktptr[j] ), 0, sizeof ( struct bucket2_r1 ) );
			memcpy ( & ( ( ( struct bucket2_r1 * ) * ( bktptr[j] ) )->kmer_t2.kmer ), & ( seq[j].kmer ), Kmer_arr_sz * sizeof ( uint64_t ) );
			( ( struct bucket2_r1 * ) * ( bktptr[j] ) )->kmer_info.cov1 = 0;
			// the cvg is useless in round 1
			pthread_spin_unlock ( &locks[hash_idx[j]] );
			( *bucket_count ) ++;
		}

		j = j + h;

		if ( j >= OverlappingKmers )
			{break;}
	}
}


/*************************************************
Function:
    process_round2_threaded_d
Description:
    Processes one read in round2:
    1. Chops read to kmers
    2. Searches the selected sparse-kmers
    3. Builds the kmer-edges (the connection between sparse-kmers)
Input:
    1. read:        a read
    2. ht:      the graph hashtable
    3. locks:       the locks array for hashtable
    4. edge_cnt:        useless
    5. K_size:      kmer size
    6. gap:         the skipped distance
Output:
    None.
Return:
    None.
*************************************************/
static void process_round2_threaded_d ( struct read_t * read, struct hashtable2 * ht, pthread_spinlock_t * locks, size_t * edge_cnt, int K_size, int gap )
{
	static size_t i;
	int readLen = read->readLen;
	int OverlappingKmers = readLen - K_size + 1;

	if ( gap >= OverlappingKmers )
		{   return;}

	int Read_arr_sz = readLen / 32 + 1;
	int rem = readLen % 32;

	if ( rem == 0 )
		{Read_arr_sz--;}

#ifdef _63MER_
	int Kmer_arr_sz = 2;
	int tot_bits = Read_arr_sz * 64;
#endif
#ifdef _127MER_
	int Kmer_arr_sz = 4;
	int tot_bits = Read_arr_sz * 128;
#endif
	size_t ht_sz = ht->ht_sz;
	bool flip[500], found[500];
	size_t hash_idx[500];
	memset ( flip, 0, sizeof ( flip ) );
	memset ( found, 0, sizeof ( found ) );
	kmer_t2 seq[500], f_seq[500];
	uint64_t hv[500], temp_bits[500];
	bucket2 ** bktptr[500];
	char c_str[500];

	for ( int j = 0; j < OverlappingKmers; j++ )
	{
		get_sub_arr ( read->read_bits, read->readLen, j, K_size, seq[j].kmer );
#ifdef _63MER_

		if ( K_size <= 31 ) //fix the represent bug
		{
			( seq[j].kmer ) [1] = ( seq[j].kmer ) [0];
			( seq[j].kmer ) [0] = 0;
		}

#endif
#ifdef _127MER_ //fix the represent bug

		if ( K_size <= 31 ) //fix the represent bug
		{
			( seq[j].kmer ) [3] = ( seq[j].kmer ) [0];
			( seq[j].kmer ) [0] = 0;
		}
		else if ( K_size <= 63 )
		{
			( seq[j].kmer ) [3] = ( seq[j].kmer ) [1];
			( seq[j].kmer ) [2] = ( seq[j].kmer ) [0];
			( seq[j].kmer ) [1] = 0;
			( seq[j].kmer ) [0] = 0;
		}
		else if ( K_size <= 95 )
		{
			( seq[j].kmer ) [3] = ( seq[j].kmer ) [2];
			( seq[j].kmer ) [2] = ( seq[j].kmer ) [1];
			( seq[j].kmer ) [1] = ( seq[j].kmer ) [0];
			( seq[j].kmer ) [0] = 0;
		}

#endif
		memcpy ( &f_seq[j], &seq[j], Kmer_arr_sz * sizeof ( uint64_t ) );
		get_rev_comp_seq_arr ( ( f_seq[j].kmer ), K_size, Kmer_arr_sz );

		if ( uint64_t_cmp ( seq[j].kmer, f_seq[j].kmer, Kmer_arr_sz ) > 0 )
		{
			flip[j] = 1;
		}

		if ( flip[j] == 1 )
		{
			memcpy ( temp_bits, & ( seq[j].kmer ), Kmer_arr_sz * sizeof ( uint64_t ) );
			memcpy ( & ( seq[j].kmer ), & ( f_seq[j].kmer ), Kmer_arr_sz * sizeof ( uint64_t ) );
			memcpy ( & ( f_seq[j].kmer ), temp_bits, Kmer_arr_sz * sizeof ( uint64_t ) );
		}

		hv[j] = MurmurHash64A ( ( seq[j].kmer ), sizeof ( seq[j] ), 0 );
		hash_idx[j] = ( size_t ) ( hv[j] % ht_sz );
		bktptr[j] = & ( ht->store_pos[hash_idx[j]] );
		found[j] = look_up_in_a_list2 ( &seq[j], &bktptr[j] );
	}

	int last_found = -1;
	int cur_found = -1;
	int h = -1;

	for ( int i = 0; i < OverlappingKmers; ++i )
	{
		if ( found[i] )
		{
			pthread_spin_lock ( &locks[hash_idx[i]] );

			if ( ( * ( bktptr[i] ) )->kmer_info.cov1 < 0xffff )
			{
				( * ( bktptr[i] ) )->kmer_info.cov1++;
			}

			pthread_spin_unlock ( &locks[hash_idx[i]] );
			cur_found = i;

			if ( last_found != -1 )
			{
				if ( cur_found - last_found > gap )
				{
					fprintf ( stderr, "ERROR: cur_found - last_found > gap !\n" );
					exit ( -1 );
				}

				//add edge ...
				h = cur_found - last_found;
				uint64_t left_bits;
				get_sub_arr ( read->read_bits, read->readLen, last_found, h, &left_bits );
				pthread_spin_lock ( &locks[hash_idx[cur_found]] ); //lock  for cur_found node to add left edge

				if ( flip[cur_found] == 0 )
				{
					struct edge_node ** edge_node_p2p = & ( ( * ( bktptr[cur_found] ) )->kmer_info.left );

					while ( ( *edge_node_p2p ) != NULL )
					{
						if ( ( *edge_node_p2p )->edge == ( uint64_t ) left_bits && ( ( *edge_node_p2p )->len + 1 ) == h )
						{
							if ( ( *edge_node_p2p )->edge_cov < 0x7f )
								{ ( *edge_node_p2p )->edge_cov++;}

							break;
						}

						edge_node_p2p = & ( ( *edge_node_p2p )->nxt_edge );
					}

					if ( ( *edge_node_p2p ) == NULL )
					{
						( *edge_node_p2p ) = ( struct edge_node * ) malloc ( sizeof ( struct edge_node ) );
						( *edge_cnt ) ++;
						memset ( *edge_node_p2p, 0, sizeof ( struct edge_node ) );
						( *edge_node_p2p )->edge = ( uint64_t ) left_bits;
						( *edge_node_p2p )->edge_cov = 1;
						( *edge_node_p2p )->len = h - 1;
					}
				}
				else
				{
					left_bits = get_rev_comp_seq ( left_bits, h );
					struct edge_node ** edge_node_p2p = & ( ( * ( bktptr[cur_found] ) )->kmer_info.right );

					while ( ( *edge_node_p2p ) != NULL )
					{
						if ( ( *edge_node_p2p )->edge == ( uint64_t ) left_bits && ( ( *edge_node_p2p )->len + 1 ) == h )
						{
							if ( ( *edge_node_p2p )->edge_cov < 0x7f )
								{ ( *edge_node_p2p )->edge_cov++;}

							break;
						}

						edge_node_p2p = & ( ( *edge_node_p2p )->nxt_edge );
					}

					if ( ( *edge_node_p2p ) == NULL )
					{
						( *edge_node_p2p ) = ( struct edge_node * ) malloc ( sizeof ( struct edge_node ) );
						( *edge_cnt ) ++;
						memset ( *edge_node_p2p, 0, sizeof ( struct edge_node ) );
						( *edge_node_p2p )->edge = ( uint64_t ) left_bits;
						( *edge_node_p2p )->edge_cov = 1;
						( *edge_node_p2p )->len = h - 1;
					}
				}

				pthread_spin_unlock ( &locks[hash_idx[cur_found]] );
				uint64_t right_bits;
				get_sub_arr ( read->read_bits, read->readLen, last_found + K_size, h, &right_bits );
				pthread_spin_lock ( &locks[hash_idx[last_found]] ); //lock ...

				if ( flip[last_found] == 1 )
				{
					right_bits = get_rev_comp_seq ( right_bits, h );
					struct edge_node ** edge_node_p2p = & ( ( * ( bktptr[last_found] ) )->kmer_info.left );

					while ( ( *edge_node_p2p ) != NULL )
					{
						if ( ( *edge_node_p2p )->edge == ( uint64_t ) right_bits && ( ( *edge_node_p2p )->len + 1 ) == h )
						{
							if ( ( *edge_node_p2p )->edge_cov < 0x7f )
								{ ( *edge_node_p2p )->edge_cov++;}

							break;
						}

						edge_node_p2p = & ( ( *edge_node_p2p )->nxt_edge );
					}

					if ( ( *edge_node_p2p ) == NULL )
					{
						( *edge_node_p2p ) = ( struct edge_node * ) malloc ( sizeof ( struct edge_node ) );
						( *edge_cnt ) ++;
						memset ( *edge_node_p2p, 0, sizeof ( struct edge_node ) );
						( *edge_node_p2p )->edge = ( uint64_t ) right_bits;
						( *edge_node_p2p )->edge_cov = 1;
						( *edge_node_p2p )->len = h - 1;
					}
				}
				else
				{
					struct edge_node ** edge_node_p2p = & ( ( * ( bktptr[last_found] ) )->kmer_info.right );

					while ( ( *edge_node_p2p ) != NULL )
					{
						if ( ( *edge_node_p2p )->edge == ( uint64_t ) right_bits && ( ( *edge_node_p2p )->len + 1 == h ) )
						{
							if ( ( *edge_node_p2p )->edge_cov < 0x7f )
								{ ( *edge_node_p2p )->edge_cov++;}

							break;
						}

						edge_node_p2p = & ( ( *edge_node_p2p )->nxt_edge );
					}

					if ( ( *edge_node_p2p ) == NULL )
					{
						( *edge_node_p2p ) = ( struct edge_node * ) malloc ( sizeof ( struct edge_node ) );
						( *edge_cnt ) ++;
						memset ( *edge_node_p2p, 0, sizeof ( struct edge_node ) );
						( *edge_node_p2p )->edge = ( uint64_t ) right_bits;
						( *edge_node_p2p )->edge_cov = 1;
						( *edge_node_p2p )->len = h - 1;
					}
				}

				pthread_spin_unlock ( &locks[hash_idx[last_found]] ); //lock ...
			}

			last_found = cur_found;
		}
	}
}


/*************************************************
Function:
    SwitchBuckets
Description:
    Switches struct bucket form round1 bucket to round2 bucket.
Input:
    1. ht:      the graph hashtable
    2. K_size:  the kmer size
Output:
    None.
Return:
    None.
*************************************************/
void SwitchBuckets ( hashtable2 * ht2, int K_size )
{
	size_t ht_sz;
	ht_sz = ht2->ht_sz;
	bucket2_r1 * store_pos_o, *store_pos_t;
	bucket2 * store_pos_n;
	bucket2 ** bktp2p;

	for ( size_t i = 0; i < ht_sz; ++i )
	{
		bktp2p = & ( ht2->store_pos[i] );
		store_pos_o = ( bucket2_r1 * ) ht2->store_pos[i];

		while ( store_pos_o != NULL )
		{
			store_pos_n = ( bucket2 * ) malloc ( sizeof ( struct bucket2 ) );
			memset ( store_pos_n, 0, sizeof ( bucket2 ) );
			store_pos_n->kmer_t2 = store_pos_o->kmer_t2;
			store_pos_n->kmer_info.cov1 = store_pos_o->kmer_info.cov1;
			*bktp2p = store_pos_n;
			bktp2p = & ( store_pos_n->nxt_bucket );
			store_pos_t = store_pos_o;
			store_pos_o = store_pos_o->nxt_bucket;
			free ( store_pos_t );
		}
	}
}


/*   c++ implementation ..

void SavingSparseKmerGraph2(hashtable2 *ht,char * outfile)
{
    ofstream  o_ht_idx,o_ht_content;

    string ht_idx_name,ht_content_name,kmer_freq;
    ht_idx_name.append(outfile).append(".ht_idx");
    ht_content_name.append(outfile).append(".ht_content");
    kmer_freq.append(outfile).append(".kmerFreq");


    map<int,int> cov_hist;

    //string ht_idx_name="HT_idx.txt",ht_content_name="HT_content";

    o_ht_idx.open(ht_idx_name.c_str(),ios_base::out|ios_base::binary);
    o_ht_content.open(ht_content_name.c_str(),ios_base::out|ios_base::binary);
    o_ht_idx<<"Hashtable size: "<<endl<<ht->ht_sz<<endl;

    bucket2 * bktptr=NULL;
    struct edge_node *edge_ptr;
    for(size_t i=0;i<ht->ht_sz;++i)
    {
        size_t list_sz=0;
        bktptr=ht->store_pos[i];
        while(bktptr!=NULL)
        {
            if(o_ht_content.write((char*) bktptr,sizeof(struct bucket2)))
            {
                edge_ptr=bktptr->kmer_info.left;
                while(edge_ptr!=NULL)
                {
                    o_ht_content.write((char*) edge_ptr,sizeof(struct edge_node));
                    edge_ptr=edge_ptr->nxt_edge;
                }
                edge_ptr=bktptr->kmer_info.right;
                while(edge_ptr!=NULL)
                {
                    o_ht_content.write((char*) edge_ptr,sizeof(struct edge_node));
                    edge_ptr=edge_ptr->nxt_edge;
                }

                int cov=bktptr->kmer_info.cov1;
                cov_hist[cov]++;
                bktptr=bktptr->nxt_bucket;
                list_sz++;
            }
            else
            {cout<<"Write error!"<<endl;}
        }
        o_ht_idx<<list_sz<<endl;
    }
    ofstream o_cov(kmer_freq.c_str());
    map<int,int >::iterator mit;
    for(mit=cov_hist.begin();mit!=cov_hist.end();++mit)
    {
        o_cov<<mit->first<<" "<<mit->second<<endl;
    }

} */

void SavingSparseKmerGraph2 ( hashtable2 * ht, char * outfile )
{
	FILE * o_ht_idx, *o_ht_content, *o_cov;
	string ht_idx_name, ht_content_name, kmer_freq;
	ht_idx_name.append ( outfile ).append ( ".ht_idx" );
	ht_content_name.append ( outfile ).append ( ".ht_content" );
	kmer_freq.append ( outfile ).append ( ".kmerFreq" );
	size_t cov_hist[256];
	memset ( cov_hist, 0, 256 * sizeof ( size_t ) );
	o_ht_idx = fopen ( ht_idx_name.c_str(), "wb" );
	o_ht_content = fopen ( ht_content_name.c_str(), "wb" );
	o_cov = fopen ( kmer_freq.c_str(), "w" );

	if ( ! ( o_ht_idx && o_ht_content && o_cov ) )
	{
		fprintf ( stderr, "ERROR: failed saving sparse kmer graph!\n" );
		return;
	}

	fprintf ( o_ht_idx, "Hashtable Size: \n" );
	fprintf ( o_ht_idx, "%llu\n", ht->ht_sz );
	bucket2 * bktptr = NULL;
	struct edge_node * edge_ptr;

	for ( size_t i = 0; i < ht->ht_sz; ++i )
	{
		size_t list_sz = 0;
		bktptr = ht->store_pos[i];

		while ( bktptr != NULL )
		{
			if ( fwrite ( ( char * ) bktptr, sizeof ( struct bucket2 ), 1, o_ht_content ) )
			{
				edge_ptr = bktptr->kmer_info.left;

				while ( edge_ptr != NULL )
				{
					fwrite ( ( char * ) edge_ptr, sizeof ( struct edge_node ), 1, o_ht_content );
					edge_ptr = edge_ptr->nxt_edge;
				}

				edge_ptr = bktptr->kmer_info.right;

				while ( edge_ptr != NULL )
				{
					fwrite ( ( char * ) edge_ptr, sizeof ( struct edge_node ), 1, o_ht_content );
					edge_ptr = edge_ptr->nxt_edge;
				}

				int cov = bktptr->kmer_info.cov1;

				if ( cov >= 255 )
				{
					cov_hist[255]++;
				}
				else
				{
					cov_hist[cov]++;
				}

				bktptr = bktptr->nxt_bucket;
				list_sz++;
			}
			else
				{cerr << "Write error!" << endl;}
		}

		fprintf ( o_ht_idx, "%llu\n", list_sz );
	}

	for ( int i = 1; i < 256; i++ )
	{
		fprintf ( o_cov, "%d\t%llu\n", i, cov_hist[i] );
	}

	fclose ( o_ht_idx );
	fclose ( o_ht_content );
	fclose ( o_cov );
}


void LoadingSparseKmerGraph2 ( hashtable2 * ht, char * outfile )
{
	string ht_idx_name, ht_content_name;
	ht_idx_name.append ( outfile ).append ( ".ht_idx" );
	ht_content_name.append ( outfile ).append ( ".ht_content" );
	ifstream in_ht_idx ( ht_idx_name.c_str(), ios_base::in | ios_base::binary ), in_ht_content ( ht_content_name.c_str(), ios_base::in | ios_base::binary );
	size_t ht_sz;
	string s;
	getline ( in_ht_idx, s );
	getline ( in_ht_idx, s );
	ht_sz = atoi ( s.c_str() ); //cerr<<ht_sz<<endl;
	Init_HT2 ( ht, ht_sz );
	struct edge_node ** edge_p2p;

	for ( int i = 0; i < ht_sz; ++i )
	{
		int list_sz;
		getline ( in_ht_idx, s );

		if ( s[s.size() - 1] == '\r' || s[s.size() - 1] == '\n' )
			{s.resize ( s.size() - 1 );}

		list_sz = atoi ( s.c_str() ); //cerr<<list_sz<<endl;
		struct bucket2 ** bktp2p = & ( ht->store_pos[i] );
		*bktp2p = NULL;

		for ( int j = 0; j < list_sz; ++j )
		{
			*bktp2p = ( struct bucket2 * ) malloc ( sizeof ( struct bucket2 ) );

			if ( in_ht_content.read ( ( char * ) ( *bktp2p ), sizeof ( struct bucket2 ) ) )
			{
				( *bktp2p )->nxt_bucket = NULL;
				( *bktp2p )->kmer_info.used = 0;
				edge_p2p = & ( ( *bktp2p )->kmer_info.left );

				while ( ( *edge_p2p ) != NULL )
				{
					( *edge_p2p ) = ( struct edge_node * ) malloc ( sizeof ( struct edge_node ) );
					in_ht_content.read ( ( char * ) ( *edge_p2p ), sizeof ( struct edge_node ) );
					edge_p2p = & ( ( *edge_p2p )->nxt_edge );
				}

				edge_p2p = & ( ( *bktp2p )->kmer_info.right );

				while ( ( *edge_p2p ) != NULL )
				{
					( *edge_p2p ) = ( struct edge_node * ) malloc ( sizeof ( struct edge_node ) );
					in_ht_content.read ( ( char * ) ( *edge_p2p ), sizeof ( struct edge_node ) );
					edge_p2p = & ( ( *edge_p2p )->nxt_edge );
				}

				bktp2p = & ( ( *bktp2p )->nxt_bucket );
			}
			else
				{cerr << "Read error!" << endl;}
		}
	}
}















































































































































































































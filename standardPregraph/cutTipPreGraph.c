/*
 * cutTipPreGraph.c
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

#include "stdinc.h"
#include "newhash.h"
#include "kmerhash.h"
#include "extfunc.h"
#include "extvab.h"

static int tip_c;  //tips counter
static long long * linearCounter; //counter for linear kmer node number


static void Mark1in1outNode ();
static void thread_mark ( KmerSet * set, unsigned char thrdID );

/*
static void printKmer(Kmer kmer)
{
    printKmerSeq(stdout,kmer);
    printf("\n");
}
*/
static int clipTipFromNode ( kmer_t * node1, int cut_len, boolean THIN )
{
	unsigned char ret = 0, in_num, out_num, link;
	int sum, count;
	kmer_t * out_node;
	Kmer tempKmer, pre_word, word, bal_word;
	ubyte8 hash_ban;
	char ch1, ch;
	boolean smaller, found;
	int setPicker;
	unsigned int max_links, singleCvg;
	/*
	if (node1->linear || node1->deleted)
	{
	    return ret;
	}

	if (THIN && !node1->single)
	{
	    return ret;
	}
	*/
	in_num = count_branch2prev ( node1 );
	out_num = count_branch2next ( node1 );

	if ( in_num == 0 && out_num == 1 )
	{
		pre_word = node1->seq;

		for ( ch1 = 0; ch1 < 4; ch1++ )
		{
			link = get_kmer_right_cov ( *node1, ch1 );

			if ( link )
			{
				break;
			}
		}

		word = nextKmer ( pre_word, ch1 );
	}
	else if ( in_num == 1 && out_num == 0 )
	{
		pre_word = reverseComplement ( node1->seq, overlaplen );

		for ( ch1 = 0; ch1 < 4; ch1++ )
		{
			link = get_kmer_left_cov ( *node1, ch1 );

			if ( link )
			{
				break;
			}
		}

		word = nextKmer ( pre_word, int_comp ( ch1 ) );
	}
	else
	{
		return ret;
	}

	count = 1;
	bal_word = reverseComplement ( word, overlaplen );

	if ( KmerLarger ( word, bal_word ) )
	{
		tempKmer = bal_word;
		bal_word = word;
		word = tempKmer;
		smaller = 0;
	}
	else
	{
		smaller = 1;
	}

	hash_ban = hash_kmer ( word );
	setPicker = hash_ban % thrd_num;
	found = search_kmerset ( KmerSets[setPicker], word, &out_node );

	if ( !found )
	{
		fprintf ( stderr, "Kmer " );
		PrintKmer ( stderr, word );
		fprintf ( stderr, " is not found, node1 " );
		PrintKmer ( stderr, node1->seq );
		fprintf ( stderr, " .\n" );
		/*
		#ifdef MER127
		fprintf (stderr,"kmer %llx%llx%llx%llx not found, node1 %llx%llx%llx%llx\n", word.high1, word.low1, word.high2, word.low2, node1->seq.high1, node1->seq.low1, node1->seq.high2,
		    node1->seq.low2);
		#else
		fprintf (stderr,"kmer %llx%llx not found, node1 %llx%llx\n", word.high, word.low, node1->seq.high, node1->seq.low);
		#endif
		*/
		exit ( 1 );
	}

	while ( out_node->linear )
	{
		count++;

		if ( THIN && !out_node->single )
		{
			break;
		}

		if ( count > cut_len )
		{
			return ret;
		}

		if ( smaller )
		{
			pre_word = word;

			for ( ch = 0; ch < 4; ch++ )
			{
				link = get_kmer_right_cov ( *out_node, ch );

				if ( link )
				{
					break;
				}
			}

			word = nextKmer ( pre_word, ch );
			bal_word = reverseComplement ( word, overlaplen );

			if ( KmerLarger ( word, bal_word ) )
			{
				tempKmer = bal_word;
				bal_word = word;
				word = tempKmer;
				smaller = 0;
			}
			else
			{
				smaller = 1;
			}

			hash_ban = hash_kmer ( word );
			setPicker = hash_ban % thrd_num;
			found = search_kmerset ( KmerSets[setPicker], word, &out_node );

			if ( !found )
			{
				fprintf ( stderr, "Kmer " );
				PrintKmer ( stderr, word );
				fprintf ( stderr, " is not found, node1 " );
				PrintKmer ( stderr, node1->seq );
				fprintf ( stderr, " .\n" );
				fprintf ( stderr, "Pre_word " );
				PrintKmer ( stderr, pre_word );
				fprintf ( stderr, " with %d(smaller).\n", ch );
				/*
				#ifdef MER127
				fprintf (stderr,"kmer %llx%llx%llx%llx not found, node1 %llx%llx%llx%llx\n", word.high1, word.low1, word.high2, word.low2, node1->seq.high1, node1->seq.low1, node1->seq.high2,
				    node1->seq.low2);
				fprintf (stderr,"pre_word %llx%llx%llx%llx with %d(smaller)\n", pre_word.high1, pre_word.low1, pre_word.high2, pre_word.low2, ch);
				#else
				fprintf (stderr,"kmer %llx%llx not found, node1 %llx%llx\n", word.high, word.low, node1->seq.high, node1->seq.low);
				fprintf (stderr,"pre_word %llx%llx with %d(smaller)\n", pre_word.high, pre_word.low, ch);
				#endif
				*/
				exit ( 1 );
			}
		}
		else
		{
			pre_word = bal_word;

			for ( ch = 0; ch < 4; ch++ )
			{
				link = get_kmer_left_cov ( *out_node, ch );

				if ( link )
				{
					break;
				}
			}

			word = nextKmer ( pre_word, int_comp ( ch ) );
			bal_word = reverseComplement ( word, overlaplen );

			if ( KmerLarger ( word, bal_word ) )
			{
				tempKmer = bal_word;
				bal_word = word;
				word = tempKmer;
				smaller = 0;
			}
			else
			{
				smaller = 1;
			}

			hash_ban = hash_kmer ( word );
			setPicker = hash_ban % thrd_num;
			found = search_kmerset ( KmerSets[setPicker], word, &out_node );

			if ( !found )
			{
				fprintf ( stderr, "Kmer " );
				PrintKmer ( stderr, word );
				fprintf ( stderr, " is not found, node1 " );
				PrintKmer ( stderr, node1->seq );
				fprintf ( stderr, " .\n" );
				fprintf ( stderr, "Pre_word " );
				PrintKmer ( stderr, reverseComplement ( pre_word, overlaplen ) );
				fprintf ( stderr, " with %d(larger).\n", int_comp ( ch ) );
				/*
				#ifdef MER127
				fprintf (stderr,"kmer %llx%llx%llx%llx not found, node1 %llx%llx%llx%llx\n", word.high1, word.low1, word.high2, word.low2, node1->seq.high1, node1->seq.low1, node1->seq.high2,
				    node1->seq.low2);
				fprintf (stderr,"pre_word %llx%llx%llx%llx with %d(larger)\n", reverseComplement (pre_word, overlaplen).high1, reverseComplement (pre_word, overlaplen).low1,
				    reverseComplement (pre_word, overlaplen).high2, reverseComplement (pre_word, overlaplen).low2, int_comp (ch));
				#else
				fprintf (stderr,"kmer %llx%llx not found, node1 %llx%llx\n", word.high, word.low, node1->seq.high, node1->seq.low);
				fprintf (stderr,"pre_word %llx%llx with %d(larger)\n", reverseComplement (pre_word, overlaplen).high, reverseComplement (pre_word, overlaplen).low, int_comp (ch));
				#endif
				*/
				exit ( 1 );
			}
		}
	}

	if ( ( sum = count_branch2next ( out_node ) + count_branch2prev ( out_node ) ) == 1 )
	{
		tip_c++;
		node1->deleted = 1;
		out_node->deleted = 1;
		return 1;
	}
	else
	{
		ch = firstCharInKmer ( pre_word );

		if ( THIN )
		{
			tip_c++;
			node1->deleted = 1;
			dislink2prevUncertain ( out_node, ch, smaller );
			out_node->linear = 0;
			return 1;
		}

		// make sure this tip doesn't provide most links to out_node
		max_links = 0;

		for ( ch1 = 0; ch1 < 4; ch1++ )
		{
			if ( smaller )
			{
				singleCvg = get_kmer_left_cov ( *out_node, ch1 );

				if ( singleCvg > max_links )
				{
					max_links = singleCvg;
				}
			}
			else
			{
				singleCvg = get_kmer_right_cov ( *out_node, ch1 );

				if ( singleCvg > max_links )
				{
					max_links = singleCvg;
				}
			}
		}

		if ( smaller && get_kmer_left_cov ( *out_node, ch ) < max_links )
		{
			tip_c++;
			node1->deleted = 1;
			dislink2prevUncertain ( out_node, ch, smaller );

			if ( count_branch2prev ( out_node ) == 1 && count_branch2next ( out_node ) == 1 )
			{
				out_node->linear = 1;
			}

			return 1;
		}

		if ( !smaller && get_kmer_right_cov ( *out_node, int_comp ( ch ) ) < max_links )
		{
			tip_c++;
			node1->deleted = 1;
			dislink2prevUncertain ( out_node, ch, smaller );

			if ( count_branch2prev ( out_node ) == 1 && count_branch2next ( out_node ) == 1 )
			{
				out_node->linear = 1;
			}

			return 1;
		}
	}

	return 0;
}


/*************************************************
Function:
    removeSingleTips
Description:
    Removes the tips starting from a kmer whose coverage is 1,
    and marks the linear kmers again (max tips length is 2 * k-mer).
Input:
    None.
Output:
    None.
Return:
    None.
*************************************************/

void removeSingleTips ()
{
	int i, flag = 0, cut_len_tip;
	kmer_t * rs;
	KmerSet * set;
	//count_ends(hash_table);
	cut_len_tip = 2 * overlaplen;   // >= maxReadLen4all-overlaplen+1 ? 2*overlaplen : maxReadLen4all-overlaplen+1;
	//if(cut_len_tip > 100) cut_len_tip = 100;
	fprintf ( stderr, "Start to remove frequency-one-kmer tips shorter than %d.\n", cut_len_tip );
	tip_c = 0;

	for ( i = 0; i < thrd_num; i++ )
	{
		set = KmerSets[i];
		set->iter_ptr = 0;

		while ( set->iter_ptr < set->size )
		{
			if ( !is_kmer_entity_null ( set->flags, set->iter_ptr ) )
			{
				rs = set->array + set->iter_ptr;

				if ( !rs->linear && !rs->deleted && rs->single )
				{
					flag += clipTipFromNode ( rs, cut_len_tip, 1 );
				}

				//              flag += clipTipFromNode (rs, cut_len_tip, 1);
			}

			set->iter_ptr++;
		}
	}

	fprintf ( stderr, "Total %d tip(s) removed.\n", tip_c );
	Mark1in1outNode ();
}


/*************************************************
Function:
    removeSingleTips
Description:
    Removes tips and marks the linear kmers again(max tips length is 2 * k-mer).
Input:
    None.
Output:
    None.
Return:
    None.
*************************************************/
void removeMinorTips ()
{
	int i, flag = 0, cut_len_tip;
	kmer_t * rs;
	KmerSet * set;
	//count_ends(hash_table);
	//cut_len_tip = 2*overlaplen >= maxReadLen4all-overlaplen+1 ? 2*overlaplen : maxReadLen4all-overlaplen+1;
	cut_len_tip = 2 * overlaplen;
	//if(cut_len_tip > 100) cut_len_tip = 100;
	fprintf ( stderr, "Start to remove tips with minority links.\n" );
	tip_c = 0;
	flag = 1;
	int round = 1;

	while ( flag )
	{
		flag = 0;

		for ( i = 0; i < thrd_num; i++ )
		{
			set = KmerSets[i];
			set->iter_ptr = 0;

			while ( set->iter_ptr < set->size )
			{
				if ( !is_kmer_entity_null ( set->flags, set->iter_ptr ) )
				{
					rs = set->array + set->iter_ptr;

					if ( !rs->linear && !rs->deleted )
					{
						flag += clipTipFromNode ( rs, cut_len_tip, 0 );
					}

					//                  flag += clipTipFromNode (rs, cut_len_tip, 0);
				}

				set->iter_ptr++;
			}

			//          fprintf (stderr,"Remove minor tips in kmer set %d is done.\n", i);
		}

		fprintf ( stderr, "%d tip(s) removed in cycle %d.\n", flag, round++ );
	}

	/*
	for (i = 0; i < thrd_num; i++)
	{
	    set = KmerSets[i];
	    flag = 1;

	    while (flag)
	    {
	        flag = 0;
	        set->iter_ptr = 0;

	        while (set->iter_ptr < set->size)
	        {
	            if (!is_kmer_entity_null (set->flags, set->iter_ptr))
	            {
	                rs = set->array + set->iter_ptr;
	                flag += clipTipFromNode (rs, cut_len_tip, 0);
	            }

	            set->iter_ptr++;
	        }
	    }

	    fprintf (stderr,"Remove minor tips in kmer set %d is done.\n", i);
	}
	*/
	fprintf ( stderr, "Total %d tip(s) removed.\n", tip_c );
	Mark1in1outNode ();
}

static void threadRoutine ( void * para )
{
	PARAMETER * prm;
	unsigned char id;
	prm = ( PARAMETER * ) para;
	id = prm->threadID;

	//printf("%dth thread with threadID %d, hash_table %p\n",id,prm.threadID,prm.hash_table);
	while ( 1 )
	{
		if ( * ( prm->selfSignal ) == 2 )
		{
			* ( prm->selfSignal ) = 0;
			break;
		}
		else if ( * ( prm->selfSignal ) == 1 )
		{
			thread_mark ( KmerSets[id], id );
			* ( prm->selfSignal ) = 0;
		}

		usleep ( 1 );
	}
}

static void creatThrds ( pthread_t * threads, PARAMETER * paras )
{
	unsigned char i;
	int temp;

	for ( i = 0; i < thrd_num; i++ )
	{
		if ( ( temp = pthread_create ( &threads[i], NULL, ( void * ) threadRoutine, & ( paras[i] ) ) ) != 0 )
		{
			fprintf ( stderr, "Create threads failed.\n" );
			exit ( 1 );
		}
	}

	fprintf ( stderr, "%d thread(s) initialized.\n", thrd_num );
}

static void thread_mark ( KmerSet * set, unsigned char thrdID )
{
	int in_num, out_num;
	kmer_t * rs;
	set->iter_ptr = 0;

	while ( set->iter_ptr < set->size )
	{
		if ( !is_kmer_entity_null ( set->flags, set->iter_ptr ) )
		{
			rs = set->array + set->iter_ptr;

			if ( rs->deleted || rs->linear )
			{
				set->iter_ptr++;
				continue;;
			}

			in_num = count_branch2prev ( rs );
			out_num = count_branch2next ( rs );

			if ( in_num == 1 && out_num == 1 )
			{
				rs->linear = 1;
				linearCounter[thrdID]++;
			}
		}

		set->iter_ptr++;
	}

	//printf("%lld more linear\n",linearCounter[thrdID]);
}

static void thread_wait ( pthread_t * threads )
{
	int i;

	for ( i = 0; i < thrd_num; i++ )
		if ( threads[i] != 0 )
		{
			pthread_join ( threads[i], NULL );
		}
}

static void sendWorkSignal ( unsigned char SIG, unsigned char * thrdSignals )
{
	int t;

	for ( t = 0; t < thrd_num; t++ )
	{
		thrdSignals[t + 1] = SIG;
	}

	while ( 1 )
	{
		usleep ( 10 );

		for ( t = 0; t < thrd_num; t++ )
			if ( thrdSignals[t + 1] )
			{
				break;
			}

		if ( t == thrd_num )
		{
			break;
		}
	}
}

static void Mark1in1outNode ()
{
	int i;
	long long counter = 0;
	pthread_t threads[thrd_num];
	unsigned char thrdSignal[thrd_num + 1];
	PARAMETER paras[thrd_num];

	for ( i = 0; i < thrd_num; i++ )
	{
		thrdSignal[i + 1] = 0;
		paras[i].threadID = i;
		paras[i].mainSignal = &thrdSignal[0];
		paras[i].selfSignal = &thrdSignal[i + 1];
	}

	creatThrds ( threads, paras );
	thrdSignal[0] = 0;
	linearCounter = ( long long * ) ckalloc ( thrd_num * sizeof ( long long ) );

	for ( i = 0; i < thrd_num; i++ )
	{
		linearCounter[i] = 0;
	}

	sendWorkSignal ( 1, thrdSignal ); //mark linear nodes
	sendWorkSignal ( 2, thrdSignal ); //stop threads
	thread_wait ( threads );

	for ( i = 0; i < thrd_num; i++ )
	{
		counter += linearCounter[i];
	}

	free ( ( void * ) linearCounter );
	fprintf ( stderr, "%lld linear node(s) marked.\n", counter );
}

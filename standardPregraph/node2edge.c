/*
 * node2edge.c
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
#include "stack.h"
#include "zlib.h"

#define KMERPTBLOCKSIZE 1000

//static int loopCounter;
static int nodeCounter;                 //count the node number including (K+1)mer node
static int edge_length_limit = 100000;  //static 'edge_seq'  nax length
static int edge_c, edgeCounter;         // current edge count number for both strand , all edge number
static preEDGE temp_edge;              // for temp use in merge_V2()
static char edge_seq[100000];       //use this static 'edge_seq ' as an temp seq in merge_V2() for speed ..

static void make_edge ( gzFile * fp );
static void merge_linearV2 ( char bal_edge, STACK * nStack, int count, gzFile * fp );
static int check_iden_kmerList ( STACK * stack1, STACK * stack2 );

//for stack
static STACK * nodeStack;           //the stack for storing linear nodes
static STACK * bal_nodeStack;       // the stack for storing the reverse complemental nodes ..


/*************************************************
Function:
    kmer2edges
Description:
    Builds edges by combining  linear kmers.
Input:
    1.outfile:      prefix of output file
Output:
    None.
Return:
    None.
*************************************************/
void kmer2edges ( char * outfile )
{
	gzFile * fp;
	char temp[256];
	sprintf ( temp, "%s.edge.gz", outfile );
	fp = gzopen ( temp, "w" );
	make_edge ( fp );
	gzclose ( fp );
	num_ed = edge_c;
}

/*************************************************
Function:
    stringBeads
Description:
    Puts the linear nodes into a stack (nodeStack).
Input:
    1. firstBead:       the first node
    2. nextch:      the next base
    3. node_c:      the node number
Output:
    None.
Return:
    None.
*************************************************/
static void stringBeads ( KMER_PT * firstBead, char nextch, int * node_c )
{
	boolean smaller, found;
	Kmer tempKmer, bal_word;
	Kmer word = firstBead->kmer;
	ubyte8 hash_ban;
	kmer_t * outgoing_node;
	int nodeCounter = 1, setPicker;
	char ch;
	unsigned char flag;
	KMER_PT * temp_pt, *prev_pt = firstBead;
	word = prev_pt->kmer;
	nodeCounter = 1;
	word = nextKmer ( word, nextch );
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
	found = search_kmerset ( KmerSets[setPicker], word, &outgoing_node );

	while ( found && ( outgoing_node->linear ) ) // for every node in this line
	{
		nodeCounter++;
		temp_pt = ( KMER_PT * ) stackPush ( nodeStack );
		temp_pt->node = outgoing_node;
		temp_pt->isSmaller = smaller;

		if ( smaller )
		{
			temp_pt->kmer = word;
		}
		else
		{
			temp_pt->kmer = bal_word;
		}

		prev_pt = temp_pt;

		if ( smaller )
		{
			for ( ch = 0; ch < 4; ch++ )
			{
				flag = get_kmer_right_cov ( *outgoing_node, ch );

				if ( flag )
				{
					break;
				}
			}

			word = nextKmer ( prev_pt->kmer, ch );
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
			found = search_kmerset ( KmerSets[setPicker], word, &outgoing_node );
		}
		else
		{
			for ( ch = 0; ch < 4; ch++ )
			{
				flag = get_kmer_left_cov ( *outgoing_node, ch );

				if ( flag )
				{
					break;
				}
			}

			word = nextKmer ( prev_pt->kmer, int_comp ( ch ) );
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
			found = search_kmerset ( KmerSets[setPicker], word, &outgoing_node );
		}
	}

	if ( outgoing_node ) //this is always true
	{
		nodeCounter++;
		temp_pt = ( KMER_PT * ) stackPush ( nodeStack );
		temp_pt->node = outgoing_node;
		temp_pt->isSmaller = smaller;

		if ( smaller )
		{
			temp_pt->kmer = word;
		}
		else
		{
			temp_pt->kmer = bal_word;
		}
	}

	*node_c = nodeCounter;
}

/*************************************************
Function:
    startEdgeFromNode
Description:
    Constructs edges from a branched node.
    For every branch (left , right):
    1. Puts the linear node into a stack.
    2. Checks whether the edge is plalindrome.
    3. Builds an edge by merging the linear nodes.
Input:
    1. node1:           the start node
    2. fp:      output file
Output:
    None.
Return:
    0.
*************************************************/
static int startEdgeFromNode ( kmer_t * node1, gzFile * fp )
{
	int node_c, palindrome;
	unsigned char flag;
	KMER_PT * ite_pt, *temp_pt;
	Kmer word1, bal_word1;
	char ch1;
	/*
	if (node1->linear || node1->deleted)
	{
	    return 0;
	}
	//   */
	// ignore floating loop
	word1 = node1->seq;
	bal_word1 = reverseComplement ( word1, overlaplen );

	// linear structure
	for ( ch1 = 0; ch1 < 4; ch1++ ) // for every node on outgoing list
	{
		flag = get_kmer_right_cov ( *node1, ch1 );

		if ( !flag )
		{
			continue;
		}

		emptyStack ( nodeStack );
		temp_pt = ( KMER_PT * ) stackPush ( nodeStack );
		temp_pt->node = node1;
		temp_pt->isSmaller = 1;
		temp_pt->kmer = word1;
		stringBeads ( temp_pt, ch1, &node_c );

		//printf("%d nodes\n",node_c);
		if ( node_c < 2 )
		{
			fprintf ( stderr, "%d nodes in this line!!!!!!!!!!!\n", node_c );
		}
		else
		{
			//make a reverse complement node list
			stackBackup ( nodeStack );
			emptyStack ( bal_nodeStack );

			while ( ( ite_pt = ( KMER_PT * ) stackPop ( nodeStack ) ) != NULL )
			{
				temp_pt = ( KMER_PT * ) stackPush ( bal_nodeStack );
				temp_pt->kmer = reverseComplement ( ite_pt->kmer, overlaplen );
			}

			stackRecover ( nodeStack );
			palindrome = check_iden_kmerList ( nodeStack, bal_nodeStack );
			stackRecover ( nodeStack );

			if ( palindrome )
			{
				merge_linearV2 ( 0, nodeStack, node_c, fp );
			}
			else
			{
				merge_linearV2 ( 1, nodeStack, node_c, fp );
			}
		}
	}           //every possible outgoing edges

	for ( ch1 = 0; ch1 < 4; ch1++ ) // for every node on incoming list
	{
		flag = get_kmer_left_cov ( *node1, ch1 );

		if ( !flag )
		{
			continue;
		}

		emptyStack ( nodeStack );
		temp_pt = ( KMER_PT * ) stackPush ( nodeStack );
		temp_pt->node = node1;
		temp_pt->isSmaller = 0;
		temp_pt->kmer = bal_word1;
		stringBeads ( temp_pt, int_comp ( ch1 ), &node_c );

		if ( node_c < 2 )
		{
			fprintf ( stderr, "%d nodes in this line!!!!!!!!!!!\n", node_c );
		}
		else
		{
			//make a reverse complement node list
			stackBackup ( nodeStack );
			emptyStack ( bal_nodeStack );

			while ( ( ite_pt = ( KMER_PT * ) stackPop ( nodeStack ) ) != NULL )
			{
				temp_pt = ( KMER_PT * ) stackPush ( bal_nodeStack );
				temp_pt->kmer = reverseComplement ( ite_pt->kmer, overlaplen );
			}

			stackRecover ( nodeStack );
			palindrome = check_iden_kmerList ( nodeStack, bal_nodeStack );
			stackRecover ( nodeStack );

			if ( palindrome )
			{
				merge_linearV2 ( 0, nodeStack, node_c, fp );
				//printf("edge is palindrome with length %d\n",temp_edge.length);
			}
			else
			{
				merge_linearV2 ( 1, nodeStack, node_c, fp );
			}
		}
	}           //every possible incoming edges

	return 0;
}

/*************************************************
Function:
    make_edge
Description:
    Builds the edges and outputs them.
Input:
    1. fp:      output file
Output:
    None.
Return:
    None.
*************************************************/
void make_edge ( gzFile * fp )
{
	int i = 0;
	kmer_t * node1;
	KmerSet * set;
	KmerSetsPatch = ( KmerSet ** ) ckalloc ( thrd_num * sizeof ( KmerSet * ) );

	for ( i = 0; i < thrd_num; i++ )
	{
		KmerSetsPatch[i] = init_kmerset ( 1000, K_LOAD_FACTOR );
	}

	nodeStack = ( STACK * ) createStack ( KMERPTBLOCKSIZE, sizeof ( KMER_PT ) );
	bal_nodeStack = ( STACK * ) createStack ( KMERPTBLOCKSIZE, sizeof ( KMER_PT ) );
	edge_c = nodeCounter = 0;
	edgeCounter = 0;

	for ( i = 0; i < thrd_num; i++ )
	{
		set = KmerSets[i];
		set->iter_ptr = 0;

		while ( set->iter_ptr < set->size )
		{
			if ( !is_kmer_entity_null ( set->flags, set->iter_ptr ) )
			{
				node1 = set->array + set->iter_ptr;

				//              /*
				if ( !node1->linear && !node1->deleted )
				{
					startEdgeFromNode ( node1, fp );
				}

				//              */
				//              startEdgeFromNode (node1, fp);
			}

			set->iter_ptr++;
		}
	}

	fprintf ( stderr, "%d (%d) edge(s) and %d extra node(s) constructed.\n", edge_c, edgeCounter, nodeCounter );
	freeStack ( nodeStack );
	freeStack ( bal_nodeStack );
}


/*************************************************
Function:
    merge_linearV2
Description:
    1. Merges the linear nodes in a stack to build an edge.
    2. Creates (K+1)mer for edges with length=1.
Input:
    1. bal_edge:        0:palindrome 1:normal
    2. nStack:      the stack of nodes
    3. count:           the number of the nodes in the nStack
    4. fp:          output file
Output:
    None.
Return:
    None.
*************************************************/
static void merge_linearV2 ( char bal_edge, STACK * nStack, int count, gzFile * fp )
{
	int length, char_index;
	preEDGE * newedge;
	kmer_t * del_node, *longNode;
	char * tightSeq, firstCh;
	long long symbol = 0;
	int len_tSeq;
	Kmer wordplus, bal_wordplus;
	ubyte8 hash_ban;
	KMER_PT * last_np = ( KMER_PT * ) stackPop ( nStack );
	KMER_PT * second_last_np = ( KMER_PT * ) stackPop ( nStack );
	KMER_PT * first_np, *second_np = NULL;
	KMER_PT * temp;
	boolean found;
	int setPicker;
	length = count - 1;
	len_tSeq = length;

	if ( len_tSeq >= edge_length_limit )
	{
		tightSeq = ( char * ) ckalloc ( len_tSeq * sizeof ( char ) );
	}
	else
	{
		tightSeq = edge_seq;
	}

	char_index = length - 1;
	newedge = &temp_edge;
	newedge->to_node = last_np->kmer;
	newedge->length = length;
	newedge->bal_edge = bal_edge;
	tightSeq[char_index--] = lastCharInKmer ( last_np->kmer );
	firstCh = firstCharInKmer ( second_last_np->kmer );
	dislink2prevUncertain ( last_np->node, firstCh, last_np->isSmaller );
	stackRecover ( nStack );

	while ( nStack->item_c > 1 )
	{
		second_np = ( KMER_PT * ) stackPop ( nStack );
	}

	first_np = ( KMER_PT * ) stackPop ( nStack );
	//unlink first node to the second one
	dislink2nextUncertain ( first_np->node, lastCharInKmer ( second_np->kmer ), first_np->isSmaller );
	//printf("from %llx, to %llx\n",first_np->node->seq,last_np->node->seq);
	//now temp is the last node in line, out_node is the second last node in line
	newedge->from_node = first_np->kmer;

	//create a long kmer for edge with length 1
	if ( length == 1 )
	{
		nodeCounter++;
		wordplus = KmerPlus ( newedge->from_node, lastCharInKmer ( newedge->to_node ) );
		bal_wordplus = reverseComplement ( wordplus, overlaplen + 1 );
		/*
		   Kmer temp = KmerPlus(reverseComplement(newedge->to_node,overlaplen),
		   lastCharInKmer(reverseComplement(newedge->from_node,overlaplen)));
		   fprintf(stderr,"(%llx %llx) (%llx %llx) (%llx %llx)\n",
		   wordplus.high,wordplus.low,temp.high,temp.low,
		   bal_wordplus.high,bal_wordplus.low);
		 */
		edge_c++;
		edgeCounter++;

		if ( KmerSmaller ( wordplus, bal_wordplus ) )
		{
			hash_ban = hash_kmer ( wordplus );
			setPicker = hash_ban % thrd_num;
			found = put_kmerset ( KmerSetsPatch[setPicker], wordplus, 4, 4, &longNode );

			if ( found )
			{
				fprintf ( stderr, "LongNode " );
				PrintKmer ( stderr, wordplus );
				fprintf ( stderr, " already exist.\n" );
				/*
				#ifdef MER127
				fprintf (stderr,"longNode %llx %llx %llx %llx already exist\n", wordplus.high1, wordplus.low1, wordplus.high2, wordplus.low2);
				#else
				fprintf (stderr,"longNode %llx %llx already exist\n", wordplus.high, wordplus.low);
				#endif
				*/
			}

			longNode->l_links = edge_c;
			longNode->twin = ( unsigned char ) ( bal_edge + 1 );
		}
		else
		{
			hash_ban = hash_kmer ( bal_wordplus );
			setPicker = hash_ban % thrd_num;
			found = put_kmerset ( KmerSetsPatch[setPicker], bal_wordplus, 4, 4, &longNode );

			if ( found )
			{
				fprintf ( stderr, "LongNode " );
				PrintKmer ( stderr, wordplus );
				fprintf ( stderr, " already exist.\n" );
				/*
				#ifdef MER127
				fprintf (stderr,"longNode %llx %llx %llx %llx already exist\n", wordplus.high1, wordplus.low1, wordplus.high2, wordplus.low2);
				#else
				fprintf (stderr,"longNode %llx %llx already exist\n", bal_wordplus.high, bal_wordplus.low);
				#endif
				*/
			}

			longNode->l_links = edge_c + bal_edge;
			longNode->twin = ( unsigned char ) ( -bal_edge + 1 );
		}
	}
	else
	{
		edge_c++;
		edgeCounter++;
	}

	stackRecover ( nStack );
	//mark all  the internal nodes
	temp = ( KMER_PT * ) stackPop ( nStack );

	while ( nStack->item_c > 1 )
	{
		temp = ( KMER_PT * ) stackPop ( nStack );
		del_node = temp->node;
		del_node->inEdge = 1;
		symbol += get_kmer_left_covs ( *del_node );
		tightSeq[char_index--] = lastCharInKmer ( temp->kmer );
	}

	stackRecover ( nStack );
	temp = ( KMER_PT * ) stackPop ( nStack );

	while ( nStack->item_c > 1 )
	{
		temp = ( KMER_PT * ) stackPop ( nStack );
		del_node = temp->node;
		del_node->inEdge = 1;

		if ( temp->isSmaller )
		{
			del_node->l_links = edge_c;
			del_node->twin = ( unsigned char ) ( bal_edge + 1 );
		}
		else
		{
			del_node->l_links = edge_c + bal_edge;
			del_node->twin = ( unsigned char ) ( -bal_edge + 1 );
		}
	}

	newedge->seq = tightSeq;

	if ( length > 1 )
	{
		newedge->cvg = symbol / ( length - 1 ) * 10 > MaxEdgeCov ? MaxEdgeCov : symbol / ( length - 1 ) * 10;
	}
	else
	{
		newedge->cvg = 0;
	}

	output_1edge ( newedge, fp );

	if ( len_tSeq >= edge_length_limit )
	{
		free ( ( void * ) tightSeq );
	}

	edge_c += bal_edge;

	if ( edge_c % 10000000 == 0 )
	{
		fprintf ( stderr, "--- %d edge(s) built.\n", edge_c );
	}

	return;
}

/*************************************************
Function:
    check_iden_kmerList
Description:
    Checks if two statcks are equal.
Input:
    1. stack1:      the first stack
    2. stack2:      the second stack
Output:
    None.
Return:
    1 if the two statcks are equal.
*************************************************/
static int check_iden_kmerList ( STACK * stack1, STACK * stack2 )
{
	KMER_PT * ite1, *ite2;

	if ( !stack1->item_c || !stack2->item_c ) // one of them is empty
	{
		return 0;
	}

	while ( ( ite1 = ( KMER_PT * ) stackPop ( stack1 ) ) != NULL && ( ite2 = ( KMER_PT * ) stackPop ( stack2 ) ) != NULL )
	{
		if ( !KmerEqual ( ite1->kmer, ite2->kmer ) )
		{
			return 0;
		}
	}

	if ( stack1->item_c || stack2->item_c ) // one of them is not empty
	{
		return 0;
	}
	else
	{
		return 1;
	}
}

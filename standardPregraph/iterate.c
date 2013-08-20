/*
 * iterate.c
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

static Kmer * kmerBuffer;
static ubyte8 * hashBanBuffer;

static char * flagBuffer;
static int buffer_size = 100000000;//2013-5-13
long long foundcount = 0;
long long notfoundcount = 0;
long long newfoundcount = 0;
long long newnotfoundcount = 0;
long long edgeaddnumber = 0;

unsigned int ** arcBuffer;
unsigned int arcBufferCount = 0;

unsigned int ** delarcBuffer;
unsigned int delarcBufferCount = 0;

static void forward ( unsigned int index, int first );

//Fresh preGraphBasic to record the original message and the multikmer message.
void freshpreGraphBasic ( boolean iterate, int maxk, char * graph )
{
	char name[256], line[256];
	FILE * fp;
	int num_kmer;
	char ch;
	int min = 0, max = 0;
	int numed = 0;
	int maxreadlen = 0, minreadlen = 0, maxnamelen = 0;
	sprintf ( name, "%s.preGraphBasic", graph );
	fp = ckopen ( name, "r" );

	while ( fgets ( line, sizeof ( line ), fp ) != NULL )
	{
		if ( line[0] == 'V' )
		{
			sscanf ( line + 6, "%d %c %d", &num_kmer, &ch, &min );
		}
		else if ( line[0] == 'E' )
		{
			sscanf ( line + 5, "%d", &numed );
		}
		else if ( line[0] == 'M' )
		{
			sscanf ( line, "MaxReadLen %d MinReadLen %d MaxNameLen %d", &maxreadlen, &minreadlen, &maxnamelen );
		}
		else if ( line[0] == 'B' )
		{
			if ( line[7] == 'V' )
			{
				sscanf ( line, "Backup VERTEX %d %c %d", &num_kmer, &ch, &min );
			}
			else if ( line[7] == 'E' )
			{
				sscanf ( line, "Backup EDGEs %d", &numed );
			}
			else if ( line[7] == 'M' )
			{
				sscanf ( line, "Backup MaxReadLen %d MinReadLen %d MaxNameLen %d", &maxreadlen, &minreadlen, &maxnamelen );
			}
		}
	}

	fclose ( fp );
	sprintf ( name, "%s.preGraphBasic", graph );
	fp = ckopen ( name, "w" );

	if ( iterate )
	{
		fprintf ( fp, "VERTEX %d K %d\n", num_vt, maxk );
		fprintf ( fp, "\nEDGEs %d\n", num_ed );
		fprintf ( fp, "\nMaxReadLen %d MinReadLen %d MaxNameLen %d\n", maxreadlen, minreadlen, maxnamelen );
		fprintf ( fp, "\nBackup VERTEX %d K %d\n", num_kmer, min );
		fprintf ( fp, "\nBackup EDGEs %d\n", numed );
		fprintf ( fp, "\nBackup MaxReadLen %d MinReadLen %d MaxNameLen %d\n", maxreadlen, minreadlen, maxnamelen );
	}
	else
	{
		fprintf ( fp, "VERTEX %d K %d\n", num_kmer, min );
		fprintf ( fp, "\nEDGEs %d\n", numed );
		fprintf ( fp, "\nMaxReadLen %d MinReadLen %d MaxNameLen %d\n", maxreadlen, minreadlen, maxnamelen );
	}

	fclose ( fp );
}

#ifdef MER127
// kmer1 = kmer1 | kmer2
inline Kmer KmerOr ( Kmer kmer1, Kmer kmer2 )
{
	kmer1.high1 |= kmer2.high1;
	kmer1.low1 |= kmer2.low1;
	kmer1.high2 |= kmer2.high2;
	kmer1.low2 |= kmer2.low2;
	return kmer1;
}

//Add ch at the head of prev.
inline Kmer KmerPlusHead ( Kmer prev, char ch, int len )
{
	Kmer word;
	word.high1 = word.low1 = word.high2 = word.low2 = 0;

	if ( 2 * len < 64 )
	{
		word.low2 = ch & 0x3;
		word.low2 <<= ( 2 * len );
	}
	else if ( 2 * len < 128 )
	{
		word.high2 = ch & 0x3;
		word.high2 <<= ( 2 * len - 64 );
	}
	else if ( 2 * len < 192 )
	{
		word.low1 = ch & 0x3;
		word.low1 <<= ( 2 * len - 128 );
	}
	else
	{
		word.high1 = ch & 0x3;
		word.high1 <<= ( 2 * len - 192 );
	}

	word = KmerOr ( word, prev );
	return word;
}


//Add ch at the tail of prev.
inline Kmer KmerPlusTail ( Kmer prev, char ch )
{
	Kmer word = KmerLeftBitMoveBy2 ( prev );
	word.low2 += ch & 0x3;
	return word;
}

static const Kmer kmerZero = { 0, 0, 0, 0 };
#else
//kmer1 = kmer1 | kmer2
Kmer KmerOr ( Kmer kmer1, Kmer kmer2 )
{
	kmer1.high |= kmer2.high;
	kmer1.low |= kmer2.low;
	return kmer1;
}

//Add ch at the head of prev.
Kmer KmerPlusHead ( Kmer prev, char ch, int len )
{
	Kmer word;
	word.high = word.low = 0;

	if ( 2 * len < 64 )
	{
		word.low = ch & 0x3;
		word.low <<= ( 2 * len );
	}
	else
	{
		word.high = ch & 0x3;
		word.high <<= ( 2 * len - 64 );
	}

	word = KmerOr ( word, prev );
	return word;
}


//Add ch at the tail of prev.
inline Kmer KmerPlusTail ( Kmer prev, char ch )
{
	Kmer word = KmerLeftBitMoveBy2 ( prev );
	word.low += ch & 0x3;
	return word;
}

static const Kmer kmerZero = { 0, 0 };
#endif

/*************************************************
Function:
    getFromKmer
Description:
    Get front (k+1)mer of edge.
Input:
    1. index :      index in edge array
Output:
    None.
Return:
    (k+1)mer at the front.
*************************************************/
Kmer  getFromKmer ( unsigned int index )
{
	Kmer temp = kmerZero;
	temp = vt_array[edge_array[index].from_vt].kmer;
	int i = 0;
	char c;

	for ( i = 0; i < step; ++i )
	{
		c = getCharInTightString ( edge_array[index].seq, i );
		temp = KmerPlusTail ( temp, c );
	}

	return temp;
}

/*************************************************
Function:
    getToKmer
Description:
    Get last (k+1)mer of edge.
Input:
    1. index :      index in edge array
Output:
    None.
Return:
    (k+1)mer at the end.
*************************************************/
Kmer  getToKmer ( unsigned int index )
{
	Kmer temp = kmerZero;
	Kmer temp2 = kmerZero;
	int len = edge_array[index].length - overlaplen + step;
	char c;
	temp = vt_array[edge_array[index].to_vt].kmer;
	temp2 = vt_array[edge_array[index].from_vt].kmer;
	int i = 0;

	if ( len >= step )
	{
		for ( i = 0; i < step; ++i )
		{
			c = getCharInTightString ( edge_array[index].seq, len - i - 1 );
			temp = KmerPlusHead ( temp, c, i + overlaplen - step );
		}
	}
	else if ( len > 0 && len < step )
	{
		for ( i = 0; i < len; ++i )
		{
			c = getCharInTightString ( edge_array[index].seq, len - i - 1 );
			temp = KmerPlusHead ( temp, c, i + overlaplen - step );
		}

		for ( i = 0; i < ( step - len ); ++i )
		{
			c = lastCharInKmer ( KmerRightBitMove ( temp2, i << 1 ) ); //.low2 & 0x3;
			temp = KmerPlusHead ( temp, c, i + overlaplen - step + len );
		}
	}
	else
	{
		for ( i = 0; i < step; ++i )
		{
			c = lastCharInKmer ( KmerRightBitMove ( temp2, ( i - len ) << 1 ) ); //.low2 & 0x3;
			temp = KmerPlusHead ( temp, c, i + overlaplen - step );
		}
	}

	return temp;
}

//Only mark on edge.
inline void delete1Edge ( unsigned int index )
{
	edge_array[index].deleted = 1;
}

/*************************************************
Function:
    kmer2vtnew
Description:
    Updates kmer edge to (k+1)mer edge.
Input:
    1. index :      edge index in array
Output:
    None.
Return:
    index of new kmer set.
    -1 if not found.
*************************************************/
int kmer2vtnew ( unsigned int index, int from )
{
	Kmer kmer;
	Kmer bal_word;
	int vt_id;
	int found = 0;
	kmer_t2 * node = NULL;

	if ( from )
		{ kmer = getFromKmer ( index ); }
	else
		{ kmer = getToKmer ( index ); }

	bal_word = reverseComplement ( kmer, overlaplen );

	if ( KmerSmaller ( kmer, bal_word ) )
	{
		vt_id = bisearch ( &vt_arraynew[0], num_vtnew, kmer );

		if ( vt_id < 0 )
		{
			++notfoundcount;
			fprintf ( stderr, "Updating edge, small vertex " );
			printKmerSeq ( stderr, kmer );
			fprintf ( stderr, " is not found, it's twin is " );
			printKmerSeq ( stderr, reverseComplement ( kmer, overlaplen ) );
			fprintf ( stderr, " .\n" );
			found = search_kmerset2 ( KmerSetsNew, kmer, &node );

			if ( found )
				{ fprintf ( stderr, "The kmer is in kmer set but not in vt_array.\n" ); }
			else
				{ fprintf ( stderr, "The kmer is not in kmer set and vt_array.\n" ); }
		}
		else
		{
			++foundcount;
		}

		return vt_id;
	}
	else
	{
		vt_id = bisearch ( &vt_arraynew[0], num_vtnew, bal_word );

		if ( vt_id >= 0 )
		{
			vt_id += num_vtnew;
			++foundcount;
		}
		else
		{
			++notfoundcount;
			fprintf ( stderr, "Updating edge, big vertex " );
			printKmerSeq ( stderr, reverseComplement ( bal_word, overlaplen ) );
			fprintf ( stderr, " is not found, it's twin is " );
			printKmerSeq ( stderr, bal_word );
			fprintf ( stderr, " .\n" );
			found = search_kmerset2 ( KmerSetsNew, bal_word, &node );

			if ( found )
				{ fprintf ( stderr, "The kmer is in kmer set but not in vt_array.\n" ); }
			else
				{ fprintf ( stderr, "The kmer is not in kmer set and vt_array.\n" ); }
		}

		return vt_id;
	}
}

/*************************************************
Function:
    update1Edge
Description:
    Updates kmer edge to (k+1)mer edge.
Input:
    1. index :      edge index in array
Output:
    None.
Return:
    None.
*************************************************/
void update1Edge ( unsigned int index )
{
	int templength = edge_array[index].length;
	int i = 0;
	int temp_from_vt;
	int temp_to_vt;
	char * tightSeq = NULL;
	temp_from_vt = kmer2vtnew ( index, 1 );
	temp_to_vt = kmer2vtnew ( index, 0 );

	if ( temp_from_vt < 0 || temp_to_vt < 0 )
	{
		destroyEdge2 ( index );
		delete1Edge ( index );
		fprintf ( stderr, "Warning : Kmer is not found, from_vt %d, to_vt %d.\n", temp_from_vt, temp_to_vt );
		return;
	}

	edge_array[index].from_vt = temp_from_vt;
	edge_array[index].to_vt = temp_to_vt;
	edge_array[index].length -= step;

	if ( edge_array[index].length == 1 || edge_array[index].length == 0 )
		{ edge_array[index].cvg = 0; }

	tightSeq = ( char * ) ckalloc ( ( edge_array[index].length / 4 + 1 ) * sizeof ( char ) );

	for ( i = 0; i < edge_array[index].length; ++i )
		{ writeChar2tightString ( getCharInTightString ( edge_array[index].seq, i + step ), tightSeq, i ); }

	if ( edge_array[index].seq )
	{
		free ( ( void * ) edge_array[index].seq );
		edge_array[index].seq = NULL;
	}

	edge_array[index].seq = tightSeq;
	edge_array[index].rv = NULL;
	ARC * currArc = edge_array[index].arcs;
	ARC * tempArc = NULL;

	while ( currArc )
	{
		tempArc = currArc;
		currArc = currArc->next;
		edge_array[index].arcs = deleteArc ( edge_array[index].arcs, tempArc );
	}

	edge_array[index].arcs = NULL;
	edge_array[index].markers = NULL;
}

/*************************************************
Function:
    getNewHash
Description:
    Gets (k+1)mer hash set.
Input:
    None.
Output:
    None.
Return:
    None.
*************************************************/
void getNewHash()
{
	unsigned int i;
	ubyte8 hash_ban, bal_hash_ban;
	Kmer word, bal_word;
	kmer_t2 * node;
	unsigned int deletecount = 0;
	kmer_cnew = 0;
	kmerBuffer = ( Kmer * ) ckalloc ( 2 * ( num_ed + 1 ) * sizeof ( Kmer ) );
	hashBanBuffer = ( ubyte8 * ) ckalloc ( 2 * ( num_ed + 1 ) * sizeof ( ubyte8 ) );
	edge_id = ( unsigned int * ) ckalloc ( 2 * ( num_ed + 1 ) * sizeof ( unsigned int ) );
	//00 : big to 01 : big from,10 : small to, 11 : small from
	flagBuffer = ( boolean * ) ckalloc ( 2 * ( num_ed + 1 ) * sizeof ( boolean ) );

	for ( i = 1; i <= num_ed; ++i )
	{
		if ( edge_array[i].deleted )
			{ continue; }

		if ( edge_array[i].length < step ) //=
		{
			destroyEdge2 ( i );
			delete1Edge ( i );
			deletecount++;
			continue;
		}

		word = kmerZero;
		bal_word = kmerZero;

		if ( edge_array[i].length == step )
		{
			word = getFromKmer ( i );
			bal_word = reverseComplement ( word, overlaplen );

			if ( KmerSmaller ( word, bal_word ) )
			{
				hash_ban = hash_kmer ( word );
				hashBanBuffer[kmer_cnew] = hash_ban;
				kmerBuffer[kmer_cnew] = word;
				flagBuffer[kmer_cnew] = 4;
			}
			else
			{
				bal_hash_ban = hash_kmer ( bal_word );
				hashBanBuffer[kmer_cnew] = bal_hash_ban;
				kmerBuffer[kmer_cnew] = bal_word;
				flagBuffer[kmer_cnew] = 5;
			}

			edge_id[kmer_cnew] = i;
			++kmer_cnew;
			continue;
		}

		if ( edge_array[i].bal_edge == 1 )
		{
			word = getFromKmer ( i );
			bal_word = reverseComplement ( word, overlaplen );

			if ( KmerSmaller ( word, bal_word ) )
			{
				hash_ban = hash_kmer ( word );
				hashBanBuffer[kmer_cnew] = hash_ban;
				kmerBuffer[kmer_cnew] = word;
				flagBuffer[kmer_cnew] = 6;
			}
			else
			{
				bal_hash_ban = hash_kmer ( bal_word );
				hashBanBuffer[kmer_cnew] = bal_hash_ban;
				kmerBuffer[kmer_cnew] = bal_word;
				flagBuffer[kmer_cnew] = 7;
			}

			edge_id[kmer_cnew] = i;
			++kmer_cnew;
			word = kmerZero;
			bal_word = kmerZero;
			word = getToKmer ( i );
			bal_word = reverseComplement ( word, overlaplen );

			if ( KmerSmaller ( word, bal_word ) )
			{
				hash_ban = hash_kmer ( word );
				hashBanBuffer[kmer_cnew] = hash_ban;
				kmerBuffer[kmer_cnew] = word;
				flagBuffer[kmer_cnew] = 8;
			}
			else
			{
				bal_hash_ban = hash_kmer ( bal_word );
				hashBanBuffer[kmer_cnew] = bal_hash_ban;
				kmerBuffer[kmer_cnew] = bal_word;
				flagBuffer[kmer_cnew] = 9;
			}

			edge_id[kmer_cnew] = i;
			++kmer_cnew;
			continue;
		}

		word = getFromKmer ( i );
		bal_word = reverseComplement ( word, overlaplen );

		if ( KmerSmaller ( word, bal_word ) )
		{
			hash_ban = hash_kmer ( word );
			hashBanBuffer[kmer_cnew] = hash_ban;
			kmerBuffer[kmer_cnew] = word;
			flagBuffer[kmer_cnew] = 3;
		}
		else
		{
			bal_hash_ban = hash_kmer ( bal_word );
			hashBanBuffer[kmer_cnew] = bal_hash_ban;
			kmerBuffer[kmer_cnew] = bal_word;
			flagBuffer[kmer_cnew] = 1;
		}

		edge_id[kmer_cnew] = i;
		++kmer_cnew;
		word = kmerZero;
		bal_word = kmerZero;
		word = getToKmer ( i );
		bal_word = reverseComplement ( word, overlaplen );

		if ( KmerSmaller ( word, bal_word ) )
		{
			hash_ban = hash_kmer ( word );
			hashBanBuffer[kmer_cnew] = hash_ban;
			kmerBuffer[kmer_cnew] = word;
			flagBuffer[kmer_cnew] = 2;
		}
		else
		{
			bal_hash_ban = hash_kmer ( bal_word );
			hashBanBuffer[kmer_cnew] = bal_hash_ban;
			kmerBuffer[kmer_cnew] = bal_word;
			flagBuffer[kmer_cnew] = 0;
		}

		edge_id[kmer_cnew] = i;
		++kmer_cnew;
	}

	fprintf ( stderr, "%lld edge(s) deleted in length of 0.\n", deletecount );
	KmerSetsNew = init_kmerset2 ( 1024, 0.77f );

	for ( i = 0; i < kmer_cnew; ++i )
	{
		put_kmerset2 ( KmerSetsNew, kmerBuffer[i], edge_id[i], flagBuffer[i], &node );
	}

	num_vtnew = count_kmerset2 ( KmerSetsNew );
	fprintf ( stderr, "%u new kmer(s).\n", num_vtnew );
	free ( kmerBuffer );
	kmerBuffer = NULL;
	free ( hashBanBuffer );
	hashBanBuffer = NULL;
	free ( edge_id );
	edge_id = NULL;
	free ( flagBuffer );
	flagBuffer = NULL;
}

/*************************************************
Function:
    getNewVertex
Description:
    Gets (k+1)mer vertex from hash.
Input:
    None.
Output:
    None.
Return:
    None.
*************************************************/
void getNewVertex()
{
	unsigned int i;
	Kmer word, bal_word;
	vt_arraynew = ( VERTEX * ) ckalloc ( 4 * ( num_vtnew + 1 ) * sizeof ( VERTEX ) );
	num_kmer_limit = 4 * num_vtnew;
	KmerSet2 * set;
	kmer_t2 * node;
	unsigned int count = 0;
	set = KmerSetsNew;
	set->iter_ptr = 0;

	while ( set->iter_ptr < set->size )
	{
		if ( !is_kmer_entity_null ( set->flags, set->iter_ptr ) )
		{
			node = set->array + set->iter_ptr;

			if ( node )
			{
				word = node->seq;
				bal_word = reverseComplement ( word, overlaplen );

				if ( KmerSmaller ( word, bal_word ) )
					{ vt_arraynew[count].kmer = word; }
				else
					{ vt_arraynew[count].kmer = bal_word; }

				++count;
			}
		}

		set->iter_ptr ++;
	}

	num_vtnew = count;
	qsort ( &vt_arraynew[0], num_vtnew, sizeof ( vt_arraynew[0] ), cmp_vertex );

	for ( i = 0; i < num_vtnew; ++i )
	{
		bal_word = reverseComplement ( vt_arraynew[i].kmer, overlaplen );
		vt_arraynew[i + num_vtnew].kmer = bal_word;
	}
}

/*************************************************
Function:
    buildGraphHash
Description:
    1. Builds (k+1)mer graph in memory based on kmer graph.
    2. Gets (k+1)mer hash set.
    3. Gets (k+1)mer vertex.
    4. Updates kmer edges to (k+1)mer edges.
Input:
    None.
Output:
    None.
Return:
    None.
*************************************************/
void buildGraphHash()
{
	unsigned int i;
	unsigned int count = 0;
	//use from kmer & to kmer to build hash
	fprintf ( stderr, "Construct new kmer hash.\n" );
	getNewHash();
	getNewVertex();
	foundcount = 0;
	notfoundcount = 0;

	for ( i = 1; i <= num_ed; ++i )
	{
		if ( edge_array[i].deleted )
			{ continue; }

		if ( edge_array[i].length < step ) //=
		{
			destroyEdge2 ( i );
			delete1Edge ( i );
			continue;
		}

		//update twin edge together
		update1Edge ( i );
		++count;
	}

	if ( notfoundcount )
	{
		fprintf ( stderr, "There are %lld kmer(s) found.\n", foundcount );
		fprintf ( stderr, "There are %lld kmer(s) not found.\n", notfoundcount );
	}

	fprintf ( stderr, "%u edge(s) updated to %dmer edge.\n", count, overlaplen );

	if ( vt_array )
	{
		free ( ( void * ) vt_array );
		vt_array = NULL;
	}

	vt_array = vt_arraynew;
	num_vt = num_vtnew;
}

//FILE *tempF = NULL;

/*************************************************
Function:
    addArc
Description:
    1. Re-build arcs by reads.
//  2. Output selected reads if -r is set(keepReadFile is true).
Input:
    1. libfile :            the reads config file
    2. graph :      the output prefix
    3. flag :           whether the selected reads file exist
    4. last :           whether it's the last iterate
    5. maxk :           the max kmer when using multikmer
//  6. keepReadFile :   keep tmp reads file that selected for building arc
Output:
    1. *.read :     the selected reads that use for assembly
Return:
    None.
*************************************************/
void addArc ( char * libfile, char * graph, int flag, int last, int maxk ) //, boolean keepReadFile
{
	/*
	if(!flag)
	{
	    char readSeqName[256];
	    sprintf(readSeqName,"%s.read",graph);
	    tempF=fopen(readSeqName,"r");
	}
	else
	{
	    tempF = NULL;
	}
	*/
	if ( flag ) // || tempF)
	{
		//      if(tempF)
		//          fclose(tempF);
		Read2edge2 ( libfile, graph, last, maxk ); //, keepReadFile
	}
	else
		{ Read2edge ( libfile, graph, maxk ); }

	unsigned int i;

	for ( i = 1; i <= num_ed; ++i )
	{
		if ( edge_array[i].deleted )
		{
			destroyEdge2 ( i );
		}
	}

	removeDeadArcs2();
}

/*************************************************
Function:
    kmer2edge
Description:
    Extend edges.
Input:
    1. from :       whether changing from kmer
    2. index :      index of edge
    3. ch :     next char to extend
    4. backup : backup char
Output:
    None.
Return:
    Index of new extended kmer.
*************************************************/
int kmer2edge ( int from, unsigned int index, char ch, char * backup )
{
	//add 1 kmer to edge
	Kmer kmer = kmerZero, bal_word = kmerZero;
	int vt_id = 0;
	char * tightSeq = NULL;

	//  if(edge_array[index].arcs && edge_array[edge_array[index].arcs->to_ed].length <= 1)
	//      return;

	if ( from <= 2 )
	{
		kmer = vt_array[edge_array[index].from_vt].kmer;
		*backup = lastCharInKmer ( kmer ); //.low2 & 3;
		kmer = prevKmer ( kmer, ch );
	}
	else
	{
		kmer = vt_array[edge_array[index].to_vt].kmer;
		*backup = ch;
		kmer = nextKmer ( kmer, ch );
	}

	bal_word = reverseComplement ( kmer, overlaplen );

	if ( KmerSmaller ( kmer, bal_word ) )
	{
		vt_id = bisearch ( &vt_array[0], num_vt, kmer );

		if ( vt_id < 0 )
		{
			++newnotfoundcount;
			fprintf ( stderr, "When extending edge 'small vertex' is not found, edge %d kmer ", index );
			printKmerSeq ( stderr, kmer );
			fprintf ( stderr, " , it's twin " );
			printKmerSeq ( stderr, bal_word );
			fprintf ( stderr, " .\n" );
		}
		else
		{
			++newfoundcount;
		}
	}
	else
	{
		vt_id = bisearch ( &vt_array[0], num_vt, bal_word );

		if ( vt_id >= 0 )
		{
			vt_id += num_vt;
			++newfoundcount;
		}
		else
		{
			++newnotfoundcount;
			fprintf ( stderr, "When extending edge 'big vertex' is not found, edge %d kmer ", index );
			printKmerSeq ( stderr, kmer );
			fprintf ( stderr, " , it's twin " );
			printKmerSeq ( stderr, bal_word );
			fprintf ( stderr, " .\n" );
		}
	}

	if ( vt_id < 0 )
	{
		return vt_id;
	}

	char backup1 = 0;
	char ch1 = int_comp ( ch );
	Kmer kmer1 = kmerZero, bal_word1 = kmerZero;
	int vt_id1 = 0;
	char * tightSeq1 = NULL;

	if ( from <= 2 )
	{
		kmer1 = vt_array[edge_array[index + 1].to_vt].kmer;
		backup1 = ch1;
		kmer1 = nextKmer ( kmer1, ch1 );
	}
	else
	{
		kmer1 = vt_array[edge_array[index + 1].from_vt].kmer;
		backup1 = lastCharInKmer ( kmer1 ); //.low2 & 3;
		kmer1 = prevKmer ( kmer1, ch1 );
	}

	bal_word1 = reverseComplement ( kmer1, overlaplen );

	if ( KmerSmaller ( kmer1, bal_word1 ) )
	{
		vt_id1 = bisearch ( &vt_array[0], num_vt, kmer1 );

		if ( vt_id1 < 0 )
		{
			++newnotfoundcount;
			fprintf ( stderr, "When extending edge 'small vertex' is not found, edge %d kmer ", index + 1 );
			printKmerSeq ( stderr, kmer1 );
			fprintf ( stderr, " , it's twin " );
			printKmerSeq ( stderr, bal_word1 );
			fprintf ( stderr, " .\n" );
		}
		else
		{
			++newfoundcount;
		}
	}
	else
	{
		vt_id1 = bisearch ( &vt_array[0], num_vt, bal_word1 );

		if ( vt_id1 >= 0 )
		{
			vt_id1 += num_vt;
			++newfoundcount;
		}
		else
		{
			++newnotfoundcount;
			fprintf ( stderr, "When extending edge big vertex is not found, edge %d kmer ", index + 1 );
			printKmerSeq ( stderr, kmer1 );
			fprintf ( stderr, " , it's twin " );
			printKmerSeq ( stderr, bal_word1 );
			fprintf ( stderr, " .\n" );
		}
	}

	if ( vt_id1 < 0 )
	{
		return vt_id1;
	}

	int i = 0;

	if ( from <= 2 )
	{
		//small
		edge_array[index].from_vt = vt_id;
		tightSeq = ( char * ) ckalloc ( ( ( edge_array[index].length + 1 ) / 4 + 1 ) * sizeof ( char ) );
		writeChar2tightString ( *backup, tightSeq, 0 );

		for ( i = 0; i < edge_array[index].length; ++i )
			{ writeChar2tightString ( getCharInTightString ( edge_array[index].seq, i ), tightSeq, i + 1 ); }

		if ( edge_array[index].seq )
		{
			free ( ( void * ) edge_array[index].seq );
			edge_array[index].seq = NULL;
		}

		edge_array[index].seq = tightSeq;
		edge_array[index].length += 1;
		//big
		edge_array[index + 1].to_vt = vt_id1;
		tightSeq1 = ( char * ) ckalloc ( ( ( edge_array[index + 1].length + 1 ) / 4 + 1 ) * sizeof ( char ) );

		for ( i = 0; i < edge_array[index + 1].length; ++i )
			{ writeChar2tightString ( getCharInTightString ( edge_array[index + 1].seq, i ), tightSeq1, i ); }

		writeChar2tightString ( backup1, tightSeq1, i );

		if ( edge_array[index + 1].seq )
		{
			free ( ( void * ) edge_array[index + 1].seq );
			edge_array[index + 1].seq = NULL;
		}

		edge_array[index + 1].seq = tightSeq1;
		edge_array[index + 1].length += 1;
	}
	else
	{
		//small
		edge_array[index].to_vt = vt_id;
		tightSeq = ( char * ) ckalloc ( ( ( edge_array[index].length + 1 ) / 4 + 1 ) * sizeof ( char ) );

		for ( i = 0; i < edge_array[index].length; ++i )
			{ writeChar2tightString ( getCharInTightString ( edge_array[index].seq, i ), tightSeq, i ); }

		writeChar2tightString ( *backup, tightSeq, i );

		if ( edge_array[index].seq )
		{
			free ( ( void * ) edge_array[index].seq );
			edge_array[index].seq = NULL;
		}

		edge_array[index].seq = tightSeq;
		edge_array[index].length += 1;
		//big
		edge_array[index + 1].from_vt = vt_id1;
		tightSeq1 = ( char * ) ckalloc ( ( ( edge_array[index + 1].length + 1 ) / 4 + 1 ) * sizeof ( char ) );
		writeChar2tightString ( backup1, tightSeq1, 0 );

		for ( i = 0; i < edge_array[index + 1].length; ++i )
			{ writeChar2tightString ( getCharInTightString ( edge_array[index + 1].seq, i ), tightSeq1, i + 1 ); }

		if ( edge_array[index + 1].seq )
		{
			free ( ( void * ) edge_array[index + 1].seq );
			edge_array[index + 1].seq = NULL;
		}

		edge_array[index + 1].seq = tightSeq1;
		edge_array[index + 1].length += 1;
	}

	return 0;
}

//Add edge and reverse complement.
void addEdge ( unsigned int from, unsigned int to, char ch, int bal_edge, unsigned int cvg )
{
	if ( num_ed_temp + 1 > num_ed_limit )
	{
		unsigned int new_num_ed = num_ed_limit * 1.2;
		edge_array = ( EDGE * ) ckrealloc ( edge_array, ( new_num_ed + 1 ) * sizeof ( EDGE ), ( num_ed_limit + 1 ) * sizeof ( EDGE ) );
		num_ed_limit = new_num_ed;
		int j;

		for ( j = num_ed_temp + 1; j <= num_ed_limit; j++ )
		{
			edge_array[j].seq = NULL;
		}

		fprintf ( stderr, "Realloc edge array.\n" );
	}

	char * tightSeq = ( char * ) ckalloc ( sizeof ( char ) );
	writeChar2tightString ( ch, tightSeq, 0 );
	edge_array[num_ed_temp + 1].from_vt = from;
	edge_array[num_ed_temp + 1].to_vt = to;
	edge_array[num_ed_temp + 1].length = 1;
	edge_array[num_ed_temp + 1].cvg = cvg;
	edge_array[num_ed_temp + 1].bal_edge = bal_edge;
	edge_array[num_ed_temp + 1].multi = 0;
	edge_array[num_ed_temp + 1].deleted = 0;
	edge_array[num_ed_temp + 1].flag = 0;

	if ( edge_array[num_ed_temp + 1].seq )
		{ free ( edge_array[num_ed_temp + 1].seq ); }

	edge_array[num_ed_temp + 1].seq = tightSeq;
	edge_array[num_ed_temp + 1].rv = NULL;
	edge_array[num_ed_temp + 1].arcs = NULL;
	edge_array[num_ed_temp + 1].markers = NULL;
	++num_ed_temp;
}

//Check whether kmers are equal to the front kmer and last kmer of the edge.
boolean checkEqual ( unsigned int from, unsigned int to, char ch, unsigned int index, unsigned int * getIndex )
{
	if ( edge_array[index].length == 1 && ( ( edge_array[index].from_vt == from && edge_array[index].to_vt == to ) ) )
	{
		return true;
	}

	return false;
}

//Whether edge exist in set.
boolean EdgeExist ( unsigned int from, unsigned int to, char ch, kmer_t2 * node, unsigned int * index )
{
	int i = 0;
	EDGEID * temp = node->edgeId;

	while ( temp )
	{
		*index = temp->edge;

		if ( checkEqual ( from, to, ch, temp->edge, index ) )
			{ return true; }

		temp = temp->next;
	}

	return false;
}

//Update edgeId in node.
void updateNode ( kmer_t2 * node, kmer_t2 * node1 )
{
	struct edgeID * edgeid;
	edgeid = ( struct edgeID * ) malloc ( sizeof ( struct edgeID ) );
	edgeid->edge = num_ed_temp + 1;
	edgeid->flag = 0;
	edgeid->next = NULL;

	if ( node->edgeId )
		{ edgeid->next = node->edgeId; }

	node->edgeId = edgeid;
	node->count++;
	edgeid = ( struct edgeID * ) malloc ( sizeof ( struct edgeID ) );
	edgeid->edge = num_ed_temp + 2;
	edgeid->flag = 1;
	edgeid->next = NULL;

	if ( node->edgeId )
		{ edgeid->next = node->edgeId; }

	node->edgeId = edgeid;
	node->count++;
	edgeid = ( struct edgeID * ) malloc ( sizeof ( struct edgeID ) );
	edgeid->edge = num_ed_temp + 1;
	edgeid->flag = 2;
	edgeid->next = NULL;

	if ( node1->edgeId )
		{ edgeid->next = node1->edgeId; }

	node1->edgeId = edgeid;
	node1->count++;
	edgeid = ( struct edgeID * ) malloc ( sizeof ( struct edgeID ) );
	edgeid->edge = num_ed_temp + 2;
	edgeid->flag = 3;
	edgeid->next = NULL;

	if ( node1->edgeId )
		{ edgeid->next = node1->edgeId; }

	node1->edgeId = edgeid;
	node1->count++;
}

/*************************************************
Function:
    checkindegree
Description:
    1. Check whether it can be solved or not.
    2. Add edge when necessary.
Input:
    1. from :       whether it's 'from kmer'
    2. from_ed :    index of 'from edge'
    3. to_ed :  index of 'to edge'
    4. node :       node of 'last kmer' of 'from edge'
    5. node1 :  node of 'front kmer' of 'to edge'
    6. maxk :       max kmer
Output:
    None.
Return:
    None.
*************************************************/
void checkindegree ( int from, unsigned int from_ed, unsigned int to_ed, kmer_t2 * node, kmer_t2 * node1, int maxk )
{
	int arcLeft_n = 0;
	int arcRight_n = 0;
	boolean exist = false;
	char ch = lastCharInKmer ( vt_array[edge_array[to_ed].from_vt].kmer ); //.low2 & 3;
	char backup = lastCharInKmer ( KmerRightBitMove ( vt_array[edge_array[from_ed].to_vt].kmer, ( overlaplen - 1 ) << 1 ) ); //.low2 & 3;
	unsigned int index;
	ARC * originalArc = NULL;

	if ( from <= 2 )
	{
		//out->in > 1
		arcCount2 ( from_ed, &arcRight_n );

		if ( arcRight_n > 1 )
		{
			exist = EdgeExist ( edge_array[from_ed].to_vt, edge_array[to_ed].from_vt, ch, node, &index );

			if ( !exist )
			{
				updateNode ( node, node1 );
				originalArc = getArcBetween ( from_ed, to_ed );
				edgeaddnumber += 2;
				addEdge ( edge_array[from_ed].to_vt, edge_array[to_ed].from_vt, ch, 2, originalArc->multiplicity * 10 );
				//              if(overlaplen + step > maxk)
				{
					arcBuffer[0][arcBufferCount] = from_ed;
					arcBuffer[1][arcBufferCount] = num_ed_temp;
					arcBuffer[2][arcBufferCount++] = to_ed;
				}
				addEdge ( edge_array[getTwinEdge ( to_ed )].to_vt, edge_array[getTwinEdge ( from_ed )].from_vt, int_comp ( backup ), 0, originalArc->multiplicity * 10 );
			}
			else
			{
				//              if(overlaplen + step > maxk)
				{
					arcBuffer[0][arcBufferCount] = from_ed;
					arcBuffer[1][arcBufferCount] = index;
					arcBuffer[2][arcBufferCount++] = to_ed;
				}
			}
		}
	}
	else
	{
		//out->in > 1
		arcCount2 ( getTwinEdge ( to_ed ), &arcLeft_n );

		if ( arcLeft_n > 1 )
		{
			exist = EdgeExist ( edge_array[from_ed].to_vt, edge_array[to_ed].from_vt, ch, node, &index );

			if ( !exist )
			{
				updateNode ( node, node1 );
				originalArc = getArcBetween ( from_ed, to_ed );
				edgeaddnumber += 2;
				addEdge ( edge_array[from_ed].to_vt, edge_array[to_ed].from_vt, ch, 2, originalArc->multiplicity * 10 );
				//              if(overlaplen + step > maxk)
				{
					arcBuffer[0][arcBufferCount] = from_ed;
					arcBuffer[1][arcBufferCount] = num_ed_temp;
					arcBuffer[2][arcBufferCount++] = to_ed;
				}
				addEdge ( edge_array[getTwinEdge ( to_ed )].to_vt, edge_array[getTwinEdge ( from_ed )].from_vt, int_comp ( backup ), 0, originalArc->multiplicity * 10 );
			}
			else
			{
				//              if(overlaplen + step > maxk)
				{
					arcBuffer[0][arcBufferCount] = from_ed;
					arcBuffer[1][arcBufferCount] = index;
					arcBuffer[2][arcBufferCount++] = to_ed;
				}
			}
		}
	}
}

//Add arc between two edges.
static void add1Arc2 ( unsigned int from_ed, unsigned int to_ed, unsigned int weight )
{
	if ( edge_array[from_ed].to_vt != edge_array[to_ed].from_vt )
	{
		fprintf ( stderr, "Warning : Inconsistant joins between %d and %d.\n", from_ed, to_ed );
	}

	unsigned int bal_fe = getTwinEdge ( from_ed );
	unsigned int bal_te = getTwinEdge ( to_ed );

	//  fprintf(stderr, "from %u, bal %u\n", from_ed, bal_fe);
	//  fprintf(stderr, "to %u, bal %u\n", to_ed, bal_te);

	if ( from_ed > num_ed_temp || to_ed > num_ed_temp || bal_fe > num_ed_temp || bal_te > num_ed_temp )
	{
		fprintf ( stderr, "Error : Edge id is out of range.\n" );
		return;
	}

	ARC * parc, *bal_parc;
	//both arcs already exist
	parc = getArcBetween ( from_ed, to_ed );

	if ( parc )
	{
		bal_parc = parc->bal_arc;
		parc->multiplicity += weight;
		bal_parc->multiplicity += weight;
		return;
	}

	//create new arcs
	parc = allocateArc ( to_ed );
	parc->multiplicity = weight;
	parc->prev = NULL;

	if ( edge_array[from_ed].arcs )
		{ edge_array[from_ed].arcs->prev = parc; }

	parc->next = edge_array[from_ed].arcs;
	edge_array[from_ed].arcs = parc;

	// A->A'
	if ( bal_te == from_ed )
	{
		parc->bal_arc = parc;
		parc->multiplicity += weight;
		return;
	}

	bal_parc = allocateArc ( bal_fe );
	bal_parc->multiplicity = weight;
	bal_parc->prev = NULL;

	if ( edge_array[bal_te].arcs )
		{ edge_array[bal_te].arcs->prev = bal_parc; }

	bal_parc->next = edge_array[bal_te].arcs;
	edge_array[bal_te].arcs = bal_parc;
	//link them to each other
	parc->bal_arc = bal_parc;
	bal_parc->bal_arc = parc;
}

//Count step between two neighbour edges.
int countstep ( unsigned int to_vt, unsigned int from_vt )
{
	if ( to_vt == from_vt )
	{
		return 0;
	}

	Kmer to, from;
	Kmer filtertemp;
	to = vt_array[to_vt].kmer;
	from = vt_array[from_vt].kmer;
	int i = 0;

	for ( i = 0; i <= nowstep2; ++i )
	{
		filtertemp = createFilter ( overlaplen - i );

		if ( KmerEqual ( KmerRightBitMove ( from, i << 1 ), KmerAnd ( to, filtertemp ) ) )
			{ return i; }
	}

	return -1;
}

/*************************************************
Function:
    freshEdge
Description:
    Extend edges or add edges.
Input:
    1. maxk :       max kmer of multi kmer
Output:
    None.
Return:
    None.
*************************************************/
void freshEdge ( int maxk )
{
	unsigned int i = 0, j = 0;
	boolean found = 0;
	kmer_t2 * node, *node1;
	int count = 0;
	char ch = 0;
	Kmer word, bal_word;
	int from_vt_id = 0, to_vt_id = 0;
	char from_backup, to_backup;
	char * tightSeq = NULL;
	int bal_ed = 0;
	int arcLeft_n = 0, arcRight_n = 0;
	ARC * temparc = NULL;
	unsigned int tempto_ed = 0;
	fprintf ( stderr, "There are %d edge(s).\n", num_ed );
	//  if(overlaplen + step > maxk)
	{
		arcBuffer = ( unsigned int ** ) ckalloc ( sizeof ( unsigned int * ) * 3 );
		arcBuffer[0] = ( unsigned int * ) ckalloc ( sizeof ( unsigned int ) * num_ed * 3 );
		arcBuffer[1] = ( unsigned int * ) ckalloc ( sizeof ( unsigned int ) * num_ed * 3 );
		arcBuffer[2] = ( unsigned int * ) ckalloc ( sizeof ( unsigned int ) * num_ed * 3 );
	}
	int count_noextend = 0;
	num_ed_temp = num_ed;

	for ( i = 1; i <= num_ed; ++i )
	{
		if ( edge_array[i].deleted || EdSameAsTwin ( i ) )
			{ continue; }

		bal_ed = getTwinEdge ( i );
		arcCount2 ( i, &arcRight_n );
		arcCount2 ( bal_ed, &arcLeft_n );

		if ( arcLeft_n == 1 )
		{
			if ( edge_array[bal_ed].to_vt != edge_array[edge_array[bal_ed].arcs->to_ed].from_vt )
			{
				ch = lastCharInKmer ( KmerRightBitMove ( vt_array[edge_array[getTwinEdge ( edge_array[bal_ed].arcs->to_ed )].to_vt].kmer, ( overlaplen - 1 ) << 1 ) ); //.low2 & 3;
				int temp = kmer2edge ( 1, i, ch, &from_backup );

				if ( temp != 0 )
					{ count_noextend++; }
			}
		}
		else if ( arcLeft_n > 1 )
		{
			word = vt_array[edge_array[i].from_vt].kmer;
			bal_word = reverseComplement ( word, overlaplen );

			if ( KmerSmaller ( word, bal_word ) )
				{ found = search_kmerset2 ( KmerSetsNew, word, &node ); }
			else
				{ found = search_kmerset2 ( KmerSetsNew, bal_word, &node ); }

			if ( !found )
			{
				fprintf ( stderr, "When refreshing edges, 'from vertex' is not found, to_vt %d kmer ", edge_array[i].from_vt );
				printKmerSeq ( stderr, vt_array[edge_array[i].from_vt].kmer );
				fprintf ( stderr, " .\n" );
				exit ( -1 );
			}

			//          if(overlaplen < maxk)
			{
				temparc = edge_array[bal_ed].arcs;

				while ( temparc )
				{
					tempto_ed = getTwinEdge ( temparc->to_ed );
					word = vt_array[edge_array[tempto_ed].to_vt].kmer;
					bal_word = reverseComplement ( word, overlaplen );

					if ( KmerSmaller ( word, bal_word ) )
						{ found = search_kmerset2 ( KmerSetsNew, word, &node1 ); }
					else
						{ found = search_kmerset2 ( KmerSetsNew, bal_word, &node1 ); }

					if ( !found )
					{
						fprintf ( stderr, "When refreshing edges, 'to vertex' is not found, to_vt %d kmer ", edge_array[tempto_ed].to_vt );
						printKmerSeq ( stderr, vt_array[edge_array[tempto_ed].to_vt].kmer );
						fprintf ( stderr, " .\n" );
						exit ( -1 );
					}

					if ( node1 && node )
						{ checkindegree ( 1, tempto_ed, i, node1, node, maxk ); }

					temparc = temparc->next;
				}
			}
		}

		if ( arcRight_n == 1 )
		{
			if ( edge_array[i].to_vt != edge_array[edge_array[i].arcs->to_ed].from_vt )
			{
				ch = lastCharInKmer ( vt_array[edge_array[edge_array[i].arcs->to_ed].from_vt].kmer ); //.low2 & 3;
				int temp = kmer2edge ( 3, i, ch, &to_backup );

				if ( temp != 0 )
					{ count_noextend++; }
			}
		}
		else if ( arcRight_n > 1 )
		{
			word = vt_array[edge_array[i].to_vt].kmer;
			bal_word = reverseComplement ( word, overlaplen );

			if ( KmerSmaller ( word, bal_word ) )
				{ found = search_kmerset2 ( KmerSetsNew, word, &node ); }
			else
				{ found = search_kmerset2 ( KmerSetsNew, bal_word, &node ); }

			if ( !found )
			{
				fprintf ( stderr, "When refreshing edges, 'to vertex' is not found, to_vt %d kmer ", edge_array[i].to_vt );
				printKmerSeq ( stderr, vt_array[edge_array[i].to_vt].kmer );
				fprintf ( stderr, " .\n" );
				exit ( -1 );
			}

			//          if(overlaplen < maxk)
			{
				temparc = edge_array[i].arcs;

				while ( temparc )
				{
					tempto_ed = temparc->to_ed;
					word = vt_array[edge_array[tempto_ed].from_vt].kmer;
					bal_word = reverseComplement ( word, overlaplen );

					if ( KmerSmaller ( word, bal_word ) )
						{ found = search_kmerset2 ( KmerSetsNew, word, &node1 ); }
					else
						{ found = search_kmerset2 ( KmerSetsNew, bal_word, &node1 ); }

					if ( !found )
					{
						fprintf ( stderr, "When refreshing edges, 'from vertex' is not found, from_vt %d kmer ", edge_array[tempto_ed].from_vt );
						printKmerSeq ( stderr, vt_array[edge_array[tempto_ed].from_vt].kmer );
						fprintf ( stderr, " .\n" );
						exit ( -1 );
					}

					if ( node1 && node )
						{ checkindegree ( 3, i, tempto_ed, node, node1, maxk ); }

					temparc = temparc->next;
				}
			}
		}

		//two edge change at the same time
		++i;
	}

	if ( count_noextend )
		{ fprintf ( stderr, "%d edge(s) not extended.\n", count_noextend ); }

	//  if(overlaplen + step > maxk)
	{
		ARC * tempArc, *tempBalArc, *originalArc, *temp, *bal_temp;
		unsigned int from = 0;
		unsigned int mid = 0;
		unsigned int to = 0;
		int count_arcdelete = 0, count_arcadd = 0;
		int arcmulti = 0;
		int arcnotfound = 0;

		for ( i = 0; i < arcBufferCount; ++i )
		{
			from = arcBuffer[0][i];
			mid = arcBuffer[1][i];
			to = arcBuffer[2][i];

			if ( from > num_ed || mid > num_ed_temp || to > num_ed )
			{
				fprintf ( stderr, "Error : Edge id is out of range.\n" );
				exit ( -1 );
			}

			originalArc = getArcBetween ( from, to );

			if ( originalArc )
			{
				arcmulti = originalArc->multiplicity;
				count_arcdelete++;
				edge_array[from].arcs = deleteArc ( edge_array[from].arcs, originalArc );
				count_arcadd += 2;
				add1Arc2 ( from, mid, arcmulti );
				add1Arc2 ( mid, to, arcmulti );
			}
			else
			{
				originalArc = getArcBetween ( getTwinEdge ( to ), getTwinEdge ( from ) );

				if ( originalArc )
				{
					arcmulti = originalArc->multiplicity;
					count_arcdelete++;
					edge_array[getTwinEdge ( to )].arcs = deleteArc ( edge_array[getTwinEdge ( to )].arcs, originalArc );
				}
				else
				{
					++arcnotfound;
					arcmulti = 2;
				}

				count_arcadd += 2;
				add1Arc2 ( from, mid, arcmulti );
				add1Arc2 ( mid, to, arcmulti );
			}
		}

		fprintf ( stderr, "Add edges to the graph: %d arc(s) deleted, %d arc(s) added.\n", count_arcdelete, count_arcadd );

		if ( arcnotfound )
			{ fprintf ( stderr, "Warning : %d arc(s) are not found when checking.\n", arcnotfound ); }

		arcBufferCount = 0;
	}
	num_ed = num_ed_temp;
	//  if(overlaplen + step > maxk)
	{
		free ( arcBuffer[2] );
		free ( arcBuffer[1] );
		free ( arcBuffer[0] );
		free ( arcBuffer );
	}
}

//Copy edge from source to target.
void copy1Edge ( EDGE * source, EDGE * target )
{
	target->from_vt = source->from_vt;
	target->to_vt = source->to_vt;
	target->length = source->length;
	target->cvg = source->cvg;
	target->multi = source->multi;

	if ( target->seq )
	{
		free ( ( void * ) target->seq );
	}

	target->seq = source->seq;
	source->seq = NULL;
	target->arcs = source->arcs;
	source->arcs = NULL;
	target->deleted = source->deleted;
}

//Check whether two bases are equal.
int BaseEqual ( char ch1, char ch2 )
{
	if ( ch1 == ch2 )
		{ return 0; }
	else if ( ch1 > ch2 )
		{ return 1; }
	else
		{ return -1; }
}

//Check whether two edges are equal.
int EdgeEqual ( unsigned int prev, unsigned int next )
{
	int i = 0;
	int length = edge_array[prev].length;
	char ch1, ch2;
	int equal = 0;

	for ( i = 0; i < length; ++i )
	{
		ch1 = int2base ( ( int ) getCharInTightString ( edge_array[prev].seq, i ) );
		ch2 = int2base ( ( int ) getCharInTightString ( edge_array[next].seq, i ) );

		if ( ( equal = BaseEqual ( ch1, ch2 ) ) )
		{
			return equal;
		}
	}

	return 0;
}
/*************************************************
Function:
    swapedge
Description:
    Re-arrange the edges, swap smaller edges at front.
Input:
    None.
Output:
    None.
Return:
    None.
*************************************************/
void swapedge()
{
	unsigned int i;
	ARC * arc, *bal_arc, *temp_arc;
	int count_swap = 0, count_equal = 0;

	for ( i = 1; i <= num_ed; ++i )
	{
		if ( edge_array[i].deleted || EdSameAsTwin ( i ) )
			{ continue; }

		if ( EdSmallerThanTwin ( i ) )
		{
			if ( KmerLarger ( vt_array[edge_array[i].from_vt].kmer, vt_array[edge_array[i + 1].from_vt].kmer ) )
			{
				count_swap++;
				copyEdge ( i, num_ed + 1 + 1 );
				copyEdge ( i + 1, num_ed + 1 );
				copyEdge ( num_ed + 1, i );
				copyEdge ( num_ed + 1 + 1, i + 1 );
				edge_array[i].bal_edge = 2;
				edge_array[i + 1].bal_edge = 0;
				//take care of the arcs
				arc = edge_array[i].arcs;

				while ( arc )
				{
					arc->bal_arc->to_ed = i + 1;
					arc = arc->next;
				}

				arc = edge_array[i + 1].arcs;

				while ( arc )
				{
					arc->bal_arc->to_ed = i;
					arc = arc->next;
				}
			}
			else if ( KmerEqual ( vt_array[edge_array[i].from_vt].kmer, vt_array[edge_array[i + 1].from_vt].kmer ) )
			{
				int temp = EdgeEqual ( i, i + 1 );

				if ( temp == 0 )
				{
					count_equal++;
					edge_array[i].bal_edge = 1;
					delete1Edge ( i + 1 );
					//take care of the arcs
					arc = edge_array[i].arcs;

					while ( arc )
					{
						arc->bal_arc->to_ed = i;
						arc = arc->next;
					}

					bal_arc = edge_array[i + 1].arcs;
					edge_array[i + 1].arcs = NULL;

					while ( bal_arc )
					{
						temp_arc = bal_arc;
						bal_arc = bal_arc->next;

						if ( edge_array[i].arcs )
							{ edge_array[i].arcs->prev = temp_arc; }

						temp_arc->next = edge_array[i].arcs;
						edge_array[i].arcs = temp_arc;
					}
				}
				else if ( temp > 0 )
				{
					count_swap++;
					copyEdge ( i, num_ed + 1 + 1 );
					copyEdge ( i + 1, num_ed + 1 );
					copyEdge ( num_ed + 1, i );
					copyEdge ( num_ed + 1 + 1, i + 1 );
					edge_array[i].bal_edge = 2;
					edge_array[i + 1].bal_edge = 0;
					//take care of the arcs
					arc = edge_array[i].arcs;

					while ( arc )
					{
						arc->bal_arc->to_ed = i + 1;
						arc = arc->next;
					}

					arc = edge_array[i + 1].arcs;

					while ( arc )
					{
						arc->bal_arc->to_ed = i;
						arc = arc->next;
					}
				}
			}

			++i;
		}
		else
		{
			delete1Edge ( i );
			fprintf ( stderr, "Warning : Front edge %d is larger than %d.\n", i, i + 1 );
		}
	}

	fprintf ( stderr, "%d none-palindrome edge(s) swapped, %d palindrome edge(s) processed.\n", count_swap, count_equal );
}
/*************************************************
Function:
    cmp_seq
Description:
    Compare two seq.
Input:
    1. a :      seq a.
    2. b :      seq b.
Output:
    None.
Return:
    1 if a larger than b.
    -1 if a smaller than b.
    0 if a equal to b.
*************************************************/
static int cmp_seq ( const void * a, const void * b )
{
	EDGE_SUB * A, *B;
	A = ( EDGE_SUB * ) a;
	B = ( EDGE_SUB * ) b;

	if ( KmerLarger ( vt_array[A->from_vt].kmer, vt_array[B->from_vt].kmer ) )
	{
		return 1;
	}
	else if ( KmerSmaller ( vt_array[A->from_vt].kmer , vt_array[B->from_vt].kmer ) )
	{
		return -1;
	}
	else
	{
		if ( A->seq[0] > B->seq[0] )
		{
			return 1;
		}
		else if ( A->seq[0] == B->seq[0] )
		{
			int i = 0;

			for ( i = 1; i < A->length && i < B->length; i++ )
			{
				if ( getCharInTightString ( A->seq, i ) > getCharInTightString ( B->seq, i ) )
					{ return 1; }
				else if ( getCharInTightString ( A->seq, i ) < getCharInTightString ( B->seq, i ) )
					{ return -1; }
			}

			if ( i == A->length && i < B->length )
				{ return -1; }
			else if ( i < A->length && i ==  B->length )
				{ return 1; }
			else
			{
				printKmerSeq ( stderr , vt_array[A->from_vt].kmer );
				fprintf ( stderr , "\n" );
				printKmerSeq ( stderr , vt_array[B->from_vt].kmer );
				fprintf ( stderr , "\n" );

				for ( i = 0; i < A->length; i++ )
				{
					fprintf ( stderr, "%c", int2base ( ( int ) getCharInTightString ( A->seq, i ) ) );
				}

				fprintf ( stderr , "\n" );

				for ( i = 0; i < B->length; i++ )
				{
					fprintf ( stderr, "%c", int2base ( ( int ) getCharInTightString ( B->seq, i ) ) );
				}

				fprintf ( stderr , "\n" );
				fprintf ( stderr, "cmp_seq:\terr\n" );
				exit ( 0 );
				return 0;
			}
		}
		else
		{
			return -1;
		}
	}
}

//Copy edge from source to target.
static void copyOneEdge ( EDGE * target , EDGE * source )
{
	target->from_vt = source->from_vt;
	target->to_vt = source->to_vt;
	target->length = source->length;
	target->cvg = source->cvg;
	target->multi = source->multi;
	target->flag = source->flag;
	target->bal_edge = source->bal_edge;
	target->seq = source->seq;
	source->seq = NULL;
	target->arcs = source->arcs;
	source->arcs = NULL ;
	target->markers = source->markers;
	source->markers = NULL;
	target->deleted = source->deleted;
}

/*************************************************
Function:
    updateArcToEd
Description:
    Update the arcs of edges.
Input:
    1. ed_index :       the index of edge
Output:
    None.
Return:
    None.
*************************************************/
static void updateArcToEd ( unsigned int ed_index )
{
	ARC * arc = edge_array[ed_index].arcs;

	while ( arc )
	{
		arc->to_ed = index_array[arc->to_ed];
		arc = arc->next;
	}
}

/*************************************************
Function:
    sortedge
Description:
    Sort edges base on seq of edges.
Input:
    None.
Output:
    None.
Return:
    None.
*************************************************/
///*
void sortedge()
{
	unsigned int index ;
	EDGE_SUB * sort_edge;
	sort_edge = ( EDGE_SUB * ) ckalloc ( sizeof ( EDGE_SUB ) * ( num_ed + 1 ) );
	unsigned int i = 1;

	for ( index = 1 ; index <= num_ed ; index ++ )
	{
		sort_edge[i].from_vt = edge_array[index].from_vt;
		sort_edge[i].seq = edge_array[index].seq;
		sort_edge[i].to_vt = index; // record old id
		sort_edge[i].length = edge_array[index].length;
		i++;

		if ( !EdSameAsTwin ( index ) )
		{
			index++;
		}
	}

	qsort ( & ( sort_edge[1] ), i - 1, sizeof ( sort_edge[1] ), cmp_seq );
	index_array = ( unsigned int * ) ckalloc ( sizeof ( unsigned int ) * ( num_ed + 1 ) ); // used to record new id
	unsigned int new_index = 1, old_index;

	for ( index = 1; index <= i - 1; index++ )
	{
		old_index = sort_edge[index].to_vt; // old id
		sort_edge[index].seq = NULL;
		index_array[old_index] = new_index++;// old id -> new id

		if ( !EdSameAsTwin ( old_index ) )
		{
			index_array[old_index + 1] = new_index++; // old id -> new id
		}
	}

	bool * copy_array = (bool * ) ckalloc ( sizeof ( bool ) * ( num_ed + 1 ) );
	EDGE *old_edge = ( EDGE * ) ckalloc ( sizeof ( EDGE ) );
	EDGE *new_edge = ( EDGE * ) ckalloc ( sizeof ( EDGE ) );
	unsigned int next_index;
	for ( index = 1; index <= num_ed; index++ )
	{
		if(!copy_array[index])
		{
			next_index = index;
			new_index = index_array[next_index];
			if(!copy_array[next_index])// && next_index != new_index
			{
				if(copy_array[new_index])
				{
					fprintf(stderr, "Copy error: never reach here.");
				}
				copy_array[next_index] = 1;
				if(next_index != new_index)
				{
					copyOneEdge (old_edge, &(edge_array[new_index]));
					copyOneEdge ( & ( edge_array[new_index] ), & ( edge_array[next_index] ) );
				}
				updateArcToEd ( new_index );

				next_index = new_index;
				new_index = index_array[next_index];
				while(!copy_array[next_index])
				{
					if(next_index == new_index)
					{
						fprintf(stderr, "Index error: never reach here.");
					}
					copy_array[next_index] = 1;
					copyOneEdge (new_edge, &(edge_array[new_index]));
					copyOneEdge ( & ( edge_array[new_index] ), old_edge);
					updateArcToEd ( new_index );
					copyOneEdge (old_edge, new_edge);

					next_index = new_index;
					new_index = index_array[next_index];
				}
			}
		}
	}

	free (copy_array);
	free (old_edge);
	free (new_edge);
	free ( index_array );
	free ( sort_edge );
	fprintf(stderr, "%d edge(s) sorted.\n", num_ed);
}
//*/
/*
void sortedge()
{
	unsigned int index ;
	EDGE * backup_edge ;
	EDGE_SUB * sort_edge;
	sort_edge = ( EDGE_SUB * ) ckalloc ( sizeof ( EDGE_SUB ) * ( num_ed + 1 ) );
	backup_edge = ( EDGE * ) ckalloc ( sizeof ( EDGE ) * ( num_ed + 1 ) );
	unsigned int i = 1;

	for ( index = 1 ; index <= num_ed ; index ++ )
	{
		sort_edge[i].from_vt = edge_array[index].from_vt;
		sort_edge[i].seq = edge_array[index].seq;
		sort_edge[i].to_vt = index; // record old id
		sort_edge[i].length = edge_array[index].length;
		i++;
		copyOneEdge ( & ( backup_edge[index] ) , & ( edge_array[index] ) );

		if ( !EdSameAsTwin ( index ) )
		{
			index++;
			copyOneEdge ( & ( backup_edge[index] ) , & ( edge_array[index] ) );
		}
	}

	qsort ( & ( sort_edge[1] ), i - 1, sizeof ( sort_edge[1] ), cmp_seq );
	index_array = ( unsigned int * ) ckalloc ( sizeof ( unsigned int ) * ( num_ed + 1 ) ); // used to record new id
	unsigned int new_index = 1, old_index;

	for ( index = 1; index <= i - 1; index++ )
	{
		old_index = sort_edge[index].to_vt; // old id
		sort_edge[index].seq = NULL;
		index_array[old_index] = new_index++;// old id -> new id

		if ( !EdSameAsTwin ( old_index ) )
		{
			index_array[old_index + 1] = new_index++; // old id -> new id
		}
	}

	for ( index = 1; index <= num_ed; index++ )
	{
		new_index = index_array[index];
		copyOneEdge ( & ( edge_array[new_index] ), & ( backup_edge[index] ) );
		updateArcToEd ( new_index );
	}

	free ( index_array );
	free ( sort_edge );
	free ( backup_edge );
	fprintf(stderr, "%d edge(s) sorted.\n", num_ed);
}
*/
/*************************************************
Function:
    delete0Edge
Description:
    Delete edge whose leght is zero.
Input:
    None.
Output:
    None.
Return:
    None.
*************************************************/
void delete0Edge()
{
	unsigned int i = 0;
	ARC * arc_left, *arc_right;
	arcBufferCount = 0;
	arcBuffer = ( unsigned int ** ) ckalloc ( sizeof ( unsigned int * ) * 3 );
	arcBuffer[0] = ( unsigned int * ) ckalloc ( sizeof ( unsigned int ) * num_ed * 3 );
	arcBuffer[1] = ( unsigned int * ) ckalloc ( sizeof ( unsigned int ) * num_ed * 3 );
	arcBuffer[2] = ( unsigned int * ) ckalloc ( sizeof ( unsigned int ) * num_ed * 3 );

	for ( i = 1; i <= num_ed; ++i )
	{
		if ( edge_array[i].deleted || EdSameAsTwin ( i ) )
			{ continue; }

		if ( edge_array[i].length == 0 )
		{
			arc_left = edge_array[i + 1].arcs;

			while ( arc_left )
			{
				arc_right = edge_array[i].arcs;

				while ( arc_right )
				{
					arcBuffer[0][arcBufferCount] = getTwinEdge ( arc_left->to_ed );
					arcBuffer[1][arcBufferCount] = ( arc_left->multiplicity + arc_right->multiplicity + 1 ) / 2;
					arcBuffer[2][arcBufferCount++] = arc_right->to_ed;
					arc_right = arc_right->next;
				}

				arc_left = arc_left->next;
			}
		}

		++i;
	}

	unsigned int from = 0;
	unsigned int multi = 0;
	unsigned int to = 0;
	int count_edgedelete = 0, count_arcadd = 0;

	for ( i = 1; i <= num_ed; ++i )
	{
		if ( edge_array[i].deleted || EdSameAsTwin ( i ) )
			{ continue; }

		if ( edge_array[i].length == 0 )
		{
			destroyEdge2 ( i );
			count_edgedelete += 2;
		}
	}

	removeDeadArcs2();

	for ( i = 0; i < arcBufferCount; ++i )
	{
		from = arcBuffer[0][i];
		multi = arcBuffer[1][i];
		to = arcBuffer[2][i];

		if ( from == 0 || to == 0 )
		{
			fprintf ( stderr, "Error : Edge id is zero.\n" );
			continue;
		}

		if ( from > num_ed || to > num_ed )
		{
			fprintf ( stderr, "Error : Edge id is out of range.\n" );
			continue;
		}

		count_arcadd++;
		add1Arc2 ( from, to, multi );
	}

	arcBufferCount = 0;
	fprintf ( stderr, "%d edge(s) in length of 0, %d arc(s) added.\n", count_edgedelete, count_arcadd );
	free ( arcBuffer[0] );
	free ( arcBuffer[1] );
	free ( arcBuffer[2] );
	free ( arcBuffer );
}

/*************************************************
Function:
    fresh
Description:
    1. Extend edges.
    2. Delete edges whose leght is zero.
    3. Re-arrange the edges, swap smaller edge at front.
Input:
    1. maxk :       max kmer of multi kmer
Output:
    None.
Return:
    None.
*************************************************/
void fresh ( int maxk )
{
	int num = 0;
	ARC * arc_temp, *parc;
	newfoundcount = 0;
	newnotfoundcount = 0;
	edgeaddnumber = 0;
	freshEdge ( maxk );
	fprintf ( stderr, "Refresh edge: %lld edge(s) added.\n", edgeaddnumber );

	if ( newnotfoundcount )
	{
		fprintf ( stderr, "Refresh edge: %d kmer(s) found.\n", newfoundcount );
		fprintf ( stderr, "Refresh edge: %d kmer(s) not found.\n", newnotfoundcount );
	}

	if ( overlaplen + step > maxk )
	{
		delete0Edge();
	}

	//swap the smaller one forward
	swapedge();
	compactEdgeArray();
}

/*************************************************
Function:
    statistics
Description:
    Output statistics of edge array, the N50 N90 longest, etc.
Input:
    1. ed_array :       the edge array
    2. ed_num :     the number of edges
Output:
    1. statistics of edges
Return:
    None.
*************************************************/
void statistics ( EDGE * ed_array, unsigned int ed_num )
{
	unsigned int i = 0;
	unsigned int * length_array;
	int flag, count, len_c;
	long long sum = 0, N90, N50;
	int signI;
	length_array = ( unsigned int * ) ckalloc ( ed_num * sizeof ( unsigned int ) );
	//first scan for number counting
	count = len_c = 0;

	for ( i = 1; i <= ed_num; i++ )
	{
		if ( ( ed_array[i].length + overlaplen - 1 ) >= len_bar )
		{
			length_array[len_c++] = ed_array[i].length + overlaplen - 1;
		}

		if ( ed_array[i].length < 1 || ed_array[i].deleted )
		{
			continue;
		}

		count++;

		if ( EdSmallerThanTwin ( i ) )
		{
			i++;
		}
	}

	sum = 0;

	for ( signI = len_c - 1; signI >= 0; signI-- )
	{
		sum += length_array[signI];
	}

	if ( len_c > 0 )
	{
		fprintf ( stderr, "\nThere are %d contig(s) longer than %d, sum up %lld bp, with average length %lld.\n", len_c, len_bar, sum, sum / len_c );
	}

	qsort ( length_array, len_c, sizeof ( length_array[0] ), cmp_int );
	fprintf ( stderr, "The longest length is %d bp, ", length_array[len_c - 1] );
	N50 = sum * 0.5;
	N90 = sum * 0.9;
	sum = flag = 0;

	for ( signI = len_c - 1; signI >= 0; signI-- )
	{
		sum += length_array[signI];

		if ( !flag && sum >= N50 )
		{
			fprintf ( stderr, "contig N50 is %d bp, ", length_array[signI] );
			flag = 1;
		}

		if ( sum >= N90 )
		{
			fprintf ( stderr, "contig N90 is %d bp.\n", length_array[signI] );
			break;
		}
	}

	free ( ( void * ) length_array );
}

/*************************************************
Function:
    sort_arc
Description:
    Sort arcs
Input:
    1. list :       unsorted arcs list
Output:
    None.
Return:
    Sorted arcs list.
*************************************************/
ARC * sort_arc ( ARC * list )
{
	if ( !list )
		{ return list; }

//	ARC * head = ( ARC * ) malloc ( sizeof ( ARC ) );
	ARC * head = ( ARC * ) ckalloc ( sizeof ( ARC ));
	head->next = list;
	list->prev = head;
	ARC * curr = list;
	ARC * temp = list;
	ARC * temp1 = NULL;

	while ( curr )
	{
		temp = curr;

		if ( temp )
		{
			temp1 = temp->next;

			while ( temp1 )
			{
				if ( temp->to_ed > temp1->to_ed )
					{ temp = temp1; }

				temp1 = temp1->next;
			}
		}

		if ( temp && temp != curr )
		{
			if ( temp->next )
			{
				temp->prev->next = temp->next;
				temp->next->prev = temp->prev;
			}
			else
			{
				temp->prev->next = NULL;
			}

			temp->next = curr;
			temp->prev = curr->prev;
			curr->prev->next = temp;
			curr->prev = temp;
		}
		else
		{
			curr = curr->next;
		}
	}

	list = head->next;
	list->prev = NULL;
	head->next = NULL;
	free ( head );
	return list;
}

//Sort disorder arcs causing by multi thread.
void freshArc()
{
	unsigned int i;
	ARC * arc_temp, *parc;

	for ( i = 1; i <= num_ed; ++i )
	{
		if ( edge_array[i].deleted )
			{ continue; }

		edge_array[i].arcs = sort_arc ( edge_array[i].arcs );
	}
	fprintf(stderr, "Arcs sorted.\n");
}

/*************************************************
Function:
    getUnlikeArc
Description:
    Delete arcs that are chose.
Input:
    None.
Output:
    None.
Return:
    None.
*************************************************/
static void deleteUnlikeArc()
{
	unsigned int i, bal;
	ARC * arc_temp;
	int count = 0;

	for ( i = 0; i < delarcBufferCount; ++i )
	{
		arc_temp = getArcBetween ( delarcBuffer[0][i], delarcBuffer[1][i] );

		if ( arc_temp )
		{
			edge_array[delarcBuffer[0][i]].arcs = deleteArc ( edge_array[delarcBuffer[0][i]].arcs, arc_temp );
			++count;
		}

		arc_temp = getArcBetween ( getTwinEdge ( delarcBuffer[1][i] ), getTwinEdge ( delarcBuffer[0][i] ) );

		if ( arc_temp )
		{
			edge_array[getTwinEdge ( delarcBuffer[1][i] )].arcs = deleteArc ( edge_array[getTwinEdge ( delarcBuffer[1][i] )].arcs, arc_temp );
			++count;
		}
	}

	fprintf ( stderr, "%d unreliable arc(s) deleted.\n", count );
	free ( delarcBuffer[0] );
	free ( delarcBuffer[1] );
	free ( delarcBuffer );
	delarcBufferCount = 0;
}

/*************************************************
Function:
    forward
Description:
    Go forward to collect the related arcs out going from index.
Input:
    1. index :      the index of edge
    2. first :      whether it's the first one to be parsed
Output:
    None.
Return:
    None.
*************************************************/
static void forward ( unsigned int index, int first )
{
	ARC * fArc, *temp;
	fArc = edge_array[index].arcs;
	unsigned int twin = getTwinEdge ( index );
	//  if(!EdSameAsTwin(index))
	{
		if ( edge_array[index].multi != 1 )
			{ edge_array[index].multi = 2; }

		if ( edge_array[twin].multi != 1 )
			{ edge_array[twin].multi = 2; }
	}
	edge_array[index].flag = 1;
	edge_array[twin].flag = 1;

	while ( fArc )
	{
		temp = fArc;
		fArc = fArc->next;
		delarcBuffer[0][delarcBufferCount] = index;
		delarcBuffer[1][delarcBufferCount++] = temp->to_ed;

		if ( edge_array[temp->to_ed].flag )
			{ continue; }

		forward ( getTwinEdge ( temp->to_ed ), 0 );
	}
}

/*************************************************
Function:
    getUnlikeArc
Description:
    Get arcs that could be processed incorrectly.
Input:
    None.
Output:
    None.
Return:
    None.
*************************************************/
static void getUnlikeArc()
{
	unsigned int i, bal;

	for ( i = 1; i <= num_ed; ++i )
	{
		if ( edge_array[i].deleted )
			{ continue; }

		if ( EdSameAsTwin ( i ) )
		{
			edge_array[i].multi = 1;
		}
	}

	delarcBuffer = ( unsigned int ** ) ckalloc ( sizeof ( unsigned int * ) * 2 );
	delarcBuffer[0] = ( unsigned int * ) ckalloc ( sizeof ( unsigned int ) * num_ed * 3 );
	delarcBuffer[1] = ( unsigned int * ) ckalloc ( sizeof ( unsigned int ) * num_ed * 3 );
	unsigned int last = 0, curr = 0;

	for ( i = 1; i <= num_ed; ++i )
	{
		if ( edge_array[i].deleted )
			{ continue; }

		if ( edge_array[i].multi == 1 )
		{
			last = delarcBufferCount;
			forward ( i, 1 );

			if ( !EdSameAsTwin ( i ) )
			{
				forward ( getTwinEdge ( i ), 1 );
				++i;
			}

			curr = delarcBufferCount;
			unsigned int j;
			edge_array[i].flag = 0;
			edge_array[getTwinEdge ( i )].flag = 0;

			for ( j = last; j < curr; ++j )
			{
				edge_array[delarcBuffer[0][j]].flag = 0;
				edge_array[getTwinEdge ( delarcBuffer[0][j] )].flag = 0;
				edge_array[delarcBuffer[1][j]].flag = 0;
				edge_array[getTwinEdge ( delarcBuffer[1][j] )].flag = 0;
			}
		}
	}
}

/*************************************************
Function:
    Iterate
Description:
    1. Iterately gets larger kmer graph and solves some repeats.
    2. Builds (k+1)mer graph based on k mer graph.
    3. Re-builds arcs by reads.
    4. Removes errors(weak edge, low cov edge and tips).
Input:
    1. libfile :                the reads config file
    2. graph :          the output prefix
    3. maxk :           the max kmer when using multikmer
//  4. keepReadFile :       keep temp reads file that selected for building arcs
    5. M :              the strength of merging bubbles
Output:
    None.
Return:
    None.
*************************************************/
void Iterate ( char * libfile, char * graph, int maxk, int M ) //boolean keepReadFile,
{
	time_t start_t, stop_t, time_bef, time_aft, inner_start, inner_stop;
	time ( &start_t );
	unsigned int i;

	for ( i = 1; i <= num_ed; ++i )
	{
		edge_array[i].multi = 0;
	}

	int cutlen = 2 * overlaplen;
	int mink = overlaplen;
	overlaplen += step;
	nowstep2 = step;
	int flag = 0;
	statistics ( edge_array, num_ed );
	fprintf ( stderr, "\nIteration start.\n" );
	int round = 1;

	while ( overlaplen <= maxk )
	{
		unsigned int j;
		time ( &inner_start );
		WORDFILTER = createFilter ( overlaplen );
		fprintf ( stderr, "\n***************************\n" );
		fprintf ( stderr, "Iteration %d, kmer: %d\n", round++, overlaplen );
		fprintf ( stderr, "Edge number: %d\n", num_ed );
		time ( &time_bef );
		//build (k+1)mer graph
		fprintf ( stderr, "Construct %dmer graph.\n", overlaplen );
		buildGraphHash();
		time ( &time_aft );
		fprintf ( stderr, "Time spent on building hash graph: %ds.\n", ( int ) ( time_aft - time_bef ) );
		time ( &time_bef );
		//add arcs for (k+1)mer graph
		fprintf ( stderr, "\nAdd arcs to graph.\n" );
		addArc ( libfile, graph, flag, maxk - overlaplen, maxk ); //, keepReadFile
		//get arcs that could be  processed incorrectly
		getUnlikeArc();
		//delete this arcs
		deleteUnlikeArc();
		flag++;
		time ( &time_aft );
		fprintf ( stderr, "Time spent on adding arcs: %ds.\n", ( int ) ( time_aft - time_bef ) );
		time ( &time_bef );
		//sort disorder arcs causing by multi thread
		fprintf ( stderr, "Sort arcs.\n" );
		freshArc();
		time ( &time_aft );
		fprintf ( stderr, "Time spent on sorting arcs: %ds.\n", ( int ) ( time_aft - time_bef ) );

		if ( deLowEdge )
		{
			time ( &time_bef );
			fprintf ( stderr, "\nRemove weak edges and low coverage edges.\n" );
			removeWeakEdges2 ( cutlen, 1, mink );
			removeLowCovEdges2 ( cutlen, deLowEdge, mink, 0 );
			time ( &time_aft );
			fprintf ( stderr, "Time spent on removing Edges: %ds\n", ( int ) ( time_aft - time_bef ) );
		}

		if ( overlaplen + step > maxk )
		{
			time ( &time_bef );
			fprintf ( stderr, "Cut tips of the graph.\n" );
			cutTipsInGraph2 ( cutlen, 0, 0 );
			time ( &time_aft );
			fprintf ( stderr, "Time spent on cutting tips: %ds.\n", ( int ) ( time_aft - time_bef ) );
		}

		time ( &time_bef );
		fprintf ( stderr, "Refresh edges.\n" );
		//refresh to extend edge and get the edge order right
		fresh ( maxk );
		time ( &time_aft );
		fprintf ( stderr, "Time spent on refreshing edges: %ds.\n", ( int ) ( time_aft - time_bef ) );
		//free kmer set
		free_kmerset2 ( KmerSetsNew );
		overlaplen += step;
		nowstep2 += step;
		statistics ( edge_array, num_ed );
		time ( &inner_stop );
		fprintf ( stderr, "Time spent on this round: %dm.\n\n", ( int ) ( inner_stop - inner_start ) / 60 );
	}

	for ( i = 1; i <= num_ed; ++i )
	{
		edge_array[i].multi = 0;
	}

	overlaplen = maxk;
	time ( &stop_t );
	fprintf ( stderr, "Iteration finished.\n" );
	fprintf ( stderr, "Time spent on iteration: %dm.\n\n", ( int ) ( stop_t - start_t ) / 60 );
}


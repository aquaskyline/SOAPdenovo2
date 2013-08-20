/*
 * inc/sparse_kmer.h
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

#ifndef _SPARSE_KMER_H
#define _SPARSE_KMER_H

#include "seq_util.h"

inline void  initKmerFilter ( int K_size, kmer_t2 * kmer_filter )
{
#ifdef _63MER_
	( kmer_filter->kmer ) [0] = 0;
	( kmer_filter->kmer ) [1] = 1LU;
	L_shift_NC ( kmer_filter->kmer, K_size * 2, 2 );

	if ( K_size <= 31 )
	{
		( kmer_filter->kmer ) [1] -= 1; //
	}
	else
	{
		( kmer_filter->kmer ) [0] -= 1; // K_size = 32 is also ok ..
		( kmer_filter->kmer ) [1] = -1; //fff..
	}

#endif
#ifdef _127MER_
	memset ( kmer_filter->kmer, 0, sizeof ( kmer_t2 ) );
	( kmer_filter->kmer ) [3] = 1LU;
	L_shift_NC ( kmer_filter->kmer, K_size * 2, 4 );

	if ( K_size <= 31 )
	{
		( kmer_filter->kmer ) [3] -= 1;
	}
	else if ( K_size <= 63 )
	{
		( kmer_filter->kmer ) [2] -= 1;
		( kmer_filter->kmer ) [3] = -1;
	}
	else if ( K_size <= 95 )
	{
		( kmer_filter->kmer ) [1] -= 1;
		( kmer_filter->kmer ) [2] = -1;
		( kmer_filter->kmer ) [3] = -1;
	}
	else if ( K_size <= 127 )
	{
		( kmer_filter->kmer ) [0] -= 1;
		( kmer_filter->kmer ) [1] = -1;
		( kmer_filter->kmer ) [2] = -1;
		( kmer_filter->kmer ) [3] = -1;
	}

#endif
}

inline void kmerAnd ( kmer_t2 * k1, kmer_t2 * k2 ) //change k1
{
#ifdef _63MER_
	( k1->kmer ) [0] &= ( k2->kmer ) [0];
	( k1->kmer ) [1] &= ( k2->kmer ) [1];
#endif
#ifdef _127MER_
	( k1->kmer ) [0] &= ( k2->kmer ) [0];
	( k1->kmer ) [1] &= ( k2->kmer ) [1];
	( k1->kmer ) [2] &= ( k2->kmer ) [2];
	( k1->kmer ) [3] &= ( k2->kmer ) [3];
#endif
}

inline void kmerOr ( kmer_t2 * k1, kmer_t2 * k2 ) //change k1
{
#ifdef _63MER_
	( k1->kmer ) [0] |= ( k2->kmer ) [0];
	( k1->kmer ) [1] |= ( k2->kmer ) [1];
#endif
#ifdef _127MER_
	( k1->kmer ) [0] |= ( k2->kmer ) [0];
	( k1->kmer ) [1] |= ( k2->kmer ) [1];
	( k1->kmer ) [2] |= ( k2->kmer ) [2];
	( k1->kmer ) [3] |= ( k2->kmer ) [3];
#endif
}

inline void kmerMoveRight ( kmer_t2 * kmer, int base_num )
{
#ifdef _63MER_
	R_shift_NC ( kmer->kmer, base_num * 2, 2 );
#endif
#ifdef _127MER_
	R_shift_NC ( kmer->kmer, base_num * 2, 4 ); // has move 32 bug
#endif
}


inline void kmerMoveLeft ( kmer_t2 * kmer, int base_num )
{
#ifdef _63MER_
	L_shift_NC ( kmer->kmer, base_num * 2, 2 );
#endif
#ifdef _127MER_
	L_shift_NC ( kmer->kmer, base_num * 2, 4 );
#endif
}

inline int kmerCompare ( kmer_t2 * k1,  kmer_t2 * k2 )
{
#ifdef _63MER_
	return uint64_t_cmp ( k1->kmer, k2->kmer, 2 );
#endif
#ifdef _127MER_
	return uint64_t_cmp ( k1->kmer, k2->kmer, 4 );
#endif
}

inline void reverseCompKmer ( kmer_t2 * kmer , int K_size ) //result stored in *kmer ...
{
#ifdef _63MER_
	get_rev_comp_seq_arr ( kmer->kmer, K_size, 2 );
#endif
#ifdef _127MER_
	get_rev_comp_seq_arr ( kmer->kmer, K_size, 4 );
#endif
}

inline bool isSmallerKmer ( kmer_t2 * kmer, int K_size )
{
	kmer_t2 f_kmer;
	f_kmer = *kmer;
	reverseCompKmer ( &f_kmer, K_size );

	if ( kmerCompare ( kmer, &f_kmer ) < 0 )
	{
		return 1;
	}

	return 0;
}

inline void  get_kmer_from_seq ( const char * seq, int len, int K_size, int pos, kmer_t2 * kmer )
{
	if ( pos + K_size > len )
	{
		fprintf ( stderr, "ERROR: get_kmer position is invalid!\n" );
		exit ( 1 );
	}

	int start = pos, end = pos + K_size, index;
	uint64_t * arr_ptr = kmer->kmer;
	memset ( arr_ptr, 0, sizeof ( kmer_t2 ) );
	int arr_sz = sizeof ( kmer_t2 ) / sizeof ( uint64_t );
	int i = 0;
	uint64_t tmp = 0;

	for ( index = start, i = 0; index < end; index++, i++ )
	{
		switch ( seq[index] )
		{
			case 'A':
				tmp = tmp << 2;
				break;
			case 'C':
				tmp = ( tmp << 2 ) | 1;
				break;
			case 'G':
				tmp = ( tmp << 2 ) | 2;
				break;
			case 'T':
				tmp = ( tmp << 2 ) | 3;
				break;
			case 'N':
				tmp = ( tmp << 2 ) | 2; //treat N as G
				break;
			default:
				tmp = ( tmp << 2 ) | 2; // treat unknown char as G, 'S'
				fprintf ( stderr, "WARNING: get_kmer_from_seq process unknown char %c\n", seq[index] );
				//exit(1);
		}

		if ( ( i + 1 ) % 32 == 0 ) //tmp is full,  tmp has stored 32 bases
		{
			arr_ptr[i / 32] = tmp;
			tmp = 0;
		}
	}

	if ( i != K_size )
	{
		fprintf ( stderr, "ERROR: i %d is K_size \n", i );
	}

	if ( K_size <= 31 )
	{
		arr_ptr[arr_sz - 1] = tmp;
	}
	else                  //if(K_size%32 != 0){  //absolute ..because K is odd
	{
		int left_bits = ( 32 - K_size % 32 ) * 2;
		tmp = tmp << left_bits;
		arr_ptr[K_size / 32] = tmp;
		kmerMoveRight ( kmer, 32 * arr_sz - K_size );
	}

	/*

	uint64_t high=0,low=0;
	int high_start=0,high_end=0;
	int low_satrt=0,low_end=0;

	if(K_size >= 33){
	    high_start = pos;
	    high_end = pos+K_size-32;

	    low_satrt = high_end;
	    low_end = low_satrt + 32;
	}else{
	    high_start = high_end = pos;
	    low_satrt = pos;
	    low_end = pos+K_size;
	}

	//debug<<"kmer ";
	for(int i=high_start;i<high_end;++i){  //dif from soapdenovo
	    //debug<<seq[i];
	    switch(seq[i]){
	        case 'A':
	            high = high<< 2;
	            break;
	        case 'C':
	            high = (high << 2)|1;
	            break;
	        case 'G':
	            high = (high << 2)|2;
	            break;
	        case 'T':
	            high = (high << 2)|3;
	            break;
	        case 'N':
	            high = (high << 2)|2;//treat N as G as same as soapdenovo
	            //debug_build<<"N occured at "<<i<<endl;
	            break;
	        default:
	            printf("error in process unknown char %c\n",seq[i]);
	            exit(1);
	    }
	}

	//debug<<" ";

	for(int i=low_satrt;i<low_end;++i){  //dif from soapdenovo
	        //debug<<seq[i];
	        switch(seq[i]){
	            case 'A':
	                low= low<< 2;
	                break;
	            case 'C':
	                low = (low << 2)|1;
	                break;
	            case 'G':
	                low = (low << 2)|2;
	                break;
	            case 'T':
	                low = (low << 2)|3;
	                break;
	            case 'N':
	                low = (low << 2)|2;//treat N as G as same as soapdenovo
	                //debug_build<<"N occured at "<<i<<endl;
	                break;
	            default:
	                printf("error in process unknown char %c\n",seq[i]);
	                exit(1);
	        }
	    }

	kmer[0]=high;
	kmer[1]=low;

	*/
}

inline void printKmer ( const kmer_t2 * kmer, FILE * fp )
{
#ifdef _63MER_
	fprintf ( fp, "%llx %llx ,\n", ( kmer->kmer ) [0], ( kmer->kmer ) [1] );
#endif
#ifdef _127MER_
	fprintf ( fp, "%llx %llx %llx %llx,\n", ( kmer->kmer ) [0], ( kmer->kmer ) [1], ( kmer->kmer ) [2], ( kmer->kmer ) [3] );
#endif
}

inline void printKmerSeq ( kmer_t2 * kmer, int K_size, FILE * fp ) //TODO printf ATCG
{
	char str[128];
	bitsarr2str ( kmer->kmer, K_size, str, sizeof ( kmer_t2 ) / sizeof ( uint64_t ) );
	fprintf ( fp, "%s ,\n", str );
}
#endif



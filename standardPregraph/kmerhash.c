/*
 * kmerhash.c
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
#include "def.h"

#define PUBLIC_FUNC
#define PROTECTED_FUNC

#define EDGEIDBLOCKSIZE 10000

#ifdef MER127
static const kmer_t2 empty_kmer2 = {{0, 0, 0, 0}, 0, 0, 0, 0};
//Get the hash key.
static inline ubyte8 modular ( KmerSet2 * set, Kmer seq )
{
	ubyte8 temp;
	temp = ( seq.high1 % set->size ) << 32 | ( seq.low1 >> 32 & 0xffffffff );
	temp = ( temp % set->size ) << 32 | ( seq.low1 & 0xffffffff );
	temp = ( temp % set->size ) << 32 | ( seq.high2 >> 32 & 0xffffffff );
	temp = ( temp % set->size ) << 32 | ( seq.high2 & 0xffffffff );
	temp = ( temp % set->size ) << 32 | ( seq.low2 >> 32 & 0xffffffff );
	temp = ( temp % set->size ) << 32 | ( seq.low2 & 0xffffffff );
	temp = ( ubyte8 ) ( temp % set->size );
	return temp;
}
#else
static const kmer_t2 empty_kmer2 = {{0, 0}, 0, 0, 0, 0};
static inline ubyte8 modular ( KmerSet2 * set, Kmer seq )
{
	ubyte8 hc;
	__uint128_t temp;
	temp = Kmer2int128 ( seq );
	hc = temp % set->size;
	return hc;
}
#endif

/*************************************************
Function:
    update_kmer2
Description:
    Updates the id and flag of the kmer.
Input:
    1. mer:     the current kmer
    2. id:      the id of edge
    3. flag:        the type of kmer(00: big "to kmer"; 01: big "from kmer"; 10: small "to kmer"; 11: small "from kmer")
Output:
    None.
Return:
    None.
*************************************************/
PUBLIC_FUNC void update_kmer2 ( kmer_t2 * mer, int id, char flag )
{
	struct edgeID * edgeid;
	edgeid = ( struct edgeID * ) malloc ( sizeof ( struct edgeID ) );
	edgeid->edge = id;
	edgeid->flag = flag;
	edgeid->next = NULL;

	if ( mer->edgeId )
		{ edgeid->next = mer->edgeId; }

	mer->edgeId = edgeid;
	mer->count++;
}

/*************************************************
Function:
    set_new_kmer2
Description:
    Initializes the new kmer.
Input:
    1. mer:     the point of the new kmer
    2. seq:     the sequence of the kmer
    3. id:      the id of edge
    4. flag:        the type of kmer(00: big "to kmer"; 01: big "from kmer"; 10: small "to kmer"; 11: small "from kmer")
Output:
    None.
Return:
    None.
*************************************************/
PUBLIC_FUNC void set_new_kmer2 ( kmer_t2 * mer, Kmer seq, int id, char flag )
{
	*mer = empty_kmer2;
	set_kmer_seq ( *mer, seq );
	update_kmer2 ( mer, id, flag );
}

//Whether it's a prime number.
static inline int is_prime_kh ( ubyte8 num )
{
	ubyte8 i, max;

	if ( num < 4 ) { return 1; }

	if ( num % 2 == 0 ) { return 0; }

	max = ( ubyte8 ) sqrt ( ( float ) num );

	for ( i = 3; i < max; i += 2 ) { if ( num % i == 0 ) { return 0; } }

	return 1;
}

//Find next prime number.
static inline ubyte8 find_next_prime_kh ( ubyte8 num )
{
	if ( num % 2 == 0 ) { num ++; }

	while ( 1 ) { if ( is_prime_kh ( num ) ) { return num; } num += 2; }
}

/*************************************************
Function:
    init_kmerset2
Description:
    Initializes the kmerset.
Input:
    1. init_size:       the initial size of the kmerset
    2. load_factor: load factor of hash
Output:
    None.
Return:
    The initial kmerset.
*************************************************/
PUBLIC_FUNC KmerSet2 * init_kmerset2 ( ubyte8 init_size, float load_factor )
{
	KmerSet2 * set;

	if ( init_size < 3 ) { init_size = 3; }
	else { init_size = find_next_prime_kh ( init_size ); }

	set = ( KmerSet2 * ) malloc ( sizeof ( KmerSet2 ) );
	set->size   = init_size;
	set->count  = 0;
	set->max    = set->size * load_factor;

	if ( load_factor <= 0 ) { load_factor = 0.25f; }
	else if ( load_factor >= 1 ) { load_factor = 0.75f; }

	set->load_factor = load_factor;
	set->iter_ptr    = 0;
	set->array = calloc ( set->size, sizeof ( kmer_t2 ) );
	set->flags = malloc ( ( set->size + 15 ) / 16 * 4 );
	memset ( set->flags, 0x55, ( set->size + 15 ) / 16 * 4 );
	return set;
}

PROTECTED_FUNC static inline ubyte8 get_kmerset2 ( KmerSet2 * set,  Kmer seq )
{
	ubyte8 hc;
	//  hc = modular (set, seq);
#ifdef MER127
	hc = modular ( set, seq );
#else
	__uint128_t temp;
	temp = Kmer2int128 ( seq );
	hc = temp % set->size;
#endif

	while ( 1 )
	{
		if ( is_kmer_entity_null ( set->flags, hc ) )
		{
			return hc;
		}
		else
		{
			if ( KmerEqual ( get_kmer_seq ( set->array[hc] ), seq ) )
				{ return hc; }
		}

		hc ++;

		if ( hc == set->size ) { hc = 0; }
	}

	return set->size;
}

/*************************************************
Function:
    search_kmerset2
Description:
    Search kmer in kmerset.
Input:
    1. set:     the kmerset
    2. seq:     the kmer
    3. rs:      the record node
Output:
    None.
Return:
    1 if found.
*************************************************/
PUBLIC_FUNC int search_kmerset2 ( KmerSet2 * set, Kmer seq, kmer_t2 ** rs )
{
	ubyte8 hc;
	//  hc = modular (set, seq);
#ifdef MER127
	hc = modular ( set, seq );
#else
	__uint128_t temp;
	temp = Kmer2int128 ( seq );
	hc = temp % set->size;
#endif

	while ( 1 )
	{
		if ( is_kmer_entity_null ( set->flags, hc ) )
		{
			return 0;
		}
		else
		{
			if ( KmerEqual ( get_kmer_seq ( set->array[hc] ), seq ) )
			{
				*rs = set->array + hc;
				return 1;
			}
		}

		hc ++;

		if ( hc == set->size ) { hc = 0; }
	}

	return 0;
}

/*************************************************
Function:
    exists_kmerset
Description:
    Whether the kmer exists in kmerset.
Input:
    1. set:     the kmerset
    2. seq:     the kmer
Output:
    None.
Return:
    1 if it exists.
*************************************************/
PUBLIC_FUNC static inline int exists_kmerset ( KmerSet2 * set, Kmer seq )
{
	ubyte8 idx;
	idx = get_kmerset2 ( set, seq );
	return !is_kmer_entity_null ( set->flags, idx );
}

/*************************************************
Function:
    encap_kmerset2
Description:
    Enlarges the kmerset if necessary.
Input:
    1. set:     the hash of kmerset
    2. num:     the element number add to kmerset
Output:
    None.
Return:
    None.
*************************************************/
PROTECTED_FUNC static inline void encap_kmerset2 ( KmerSet2 * set, ubyte8 num )
{
	ubyte4 * flags, *f;
	ubyte8 i, n, size, hc;
	kmer_t2 key, tmp;

	if ( set->count + num <= set->max ) { return; }

	if ( initKmerSetSize != 0 )
	{
		if ( set->load_factor < 0.88 )
		{
			set->load_factor = 0.88;
			set->max    = set->size * set->load_factor;
			return;
		}
		else
		{
			fprintf ( stderr, "-- Static memory pool exploded, please define a larger value. --\nloadFactor\t%f\nsize\t%llu\ncnt\t%llu\n", set->load_factor, set->size, set->count );
			abort();
		}
	}

	n = set->size;

	do
	{
		if ( n < 0xFFFFFFFU )
			{ n <<= 1; }
		else
			{ n += 0xFFFFFFU; }

		n = find_next_prime_kh ( n );
	}
	while ( n * set->load_factor < set->count + num );

	set->array = realloc ( set->array, n * sizeof ( kmer_t2 ) );

	if ( set->array == NULL )
	{
		fprintf ( stderr, "-- Out of memory --\n" );
		abort();
	}

	flags = malloc ( ( n + 15 ) / 16 * 4 );
	memset ( flags, 0x55, ( n + 15 ) / 16 * 4 );
	size = set->size;
	set->size = n;
	set->max = n * set->load_factor;
	f = set->flags;
	set->flags = flags;
	flags = f;
	__uint128_t temp;

	for ( i = 0; i < size; i++ )
	{
		if ( !exists_kmer_entity ( flags, i ) ) { continue; }

		key = set->array[i];
		set_kmer_entity_del ( flags, i );

		while ( 1 )
		{
			hc = modular ( set, get_kmer_seq ( key ) );
#ifdef MER127
			hc = modular ( set, get_kmer_seq ( key ) );
#else
			temp = Kmer2int128 ( get_kmer_seq ( key ) );
			hc = temp % set->size;
#endif

			while ( !is_kmer_entity_null ( set->flags, hc ) ) { hc ++; if ( hc == set->size ) { hc = 0; } }

			clear_kmer_entity_null ( set->flags, hc );

			if ( hc < size && exists_kmer_entity ( flags, hc ) )
			{
				tmp = key;
				key = set->array[hc];
				set->array[hc] = tmp;
				set_kmer_entity_del ( flags, hc );
			}
			else
			{
				set->array[hc] = key;
				break;
			}
		}
	}

	free ( flags );
}

/*************************************************
Function:
    put_kmerset2
Description:
    Puts kmer into the kmerset.
Input:
    1. set:     the hash of kmerset
    2. seq:     the sequence of the kmer
    3. id:      the id of edge
    4. flag:        the type of kmer(00: big "to kmer"; 01: big "from kmer"; 10: small "to kmer"; 11: small "from kmer")
    5. kmer_p:  used to record the info of the kmer
Output:
    None.
Return:
    1 if it's successfully put kmer into kmerset.
*************************************************/
PUBLIC_FUNC int put_kmerset2 ( KmerSet2 * set, Kmer seq, int id, char flag, kmer_t2 ** kmer_p )
{
	ubyte8 hc;

	if ( set->count + 1 > set->max )
	{
		encap_kmerset2 ( set, 1 );
	}

	//  hc = modular (set, seq);
#ifdef MER127
	hc = modular ( set, seq );
#else
	__uint128_t temp;
	temp = Kmer2int128 ( seq );
	hc = temp % set->size;
#endif

	do
	{
		if ( is_kmer_entity_null ( set->flags, hc ) )
		{
			clear_kmer_entity_null ( set->flags, hc );
			set_new_kmer2 ( set->array + hc, seq, id, flag );
			set->count ++;
			*kmer_p = set->array + hc;
			return 0;
		}
		else
		{
			if ( KmerEqual ( get_kmer_seq ( set->array[hc] ), seq ) )
			{
				update_kmer2 ( set->array + hc, id, flag );
				*kmer_p = set->array + hc;
				return 1;
			}
		}

		hc ++;

		if ( hc == set->size ) { hc = 0; }
	}
	while ( 1 );

	*kmer_p = NULL;
	return 0;
}

/*************************************************
Function:
    count_kmerset2
Description:
    Returns the kmer number of the kmerset.
Input:
    None.
Output:
    None.
Return:
    The kmer number of the kmerset.
*************************************************/
PUBLIC_FUNC byte8 count_kmerset2 ( KmerSet2 * set ) { return set->count; }

PUBLIC_FUNC static inline void reset_iter_kmerset2 ( KmerSet2 * set ) { set->iter_ptr = 0; }

PUBLIC_FUNC static inline ubyte8 iter_kmerset2 ( KmerSet2 * set, kmer_t2 ** rs )
{
	while ( set->iter_ptr < set->size )
	{
		if ( !is_kmer_entity_null ( set->flags, set->iter_ptr ) )
		{
			*rs = set->array + set->iter_ptr;
			set->iter_ptr ++;
			return 1;
		}

		set->iter_ptr ++;
	}

	return 0;
}

//Free.
PUBLIC_FUNC void free_kmerset2 ( KmerSet2 * set )
{
	int i;
	struct edgeID * temp, *temp_next;
	kmer_t2 * node;
	set->iter_ptr = 0;

	while ( set->iter_ptr < set->size )
	{
		if ( !is_kmer_entity_null ( set->flags, set->iter_ptr ) )
		{
			node = set->array + set->iter_ptr;

			if ( node )
			{
				temp = node->edgeId;

				while ( temp )
				{
					temp_next = temp->next;
					free ( ( void * ) temp );
					temp = temp_next;
				}
			}
		}

		set->iter_ptr ++;
	}

	free ( set->array );
	set->array = NULL;
	free ( set->flags );
	set->flags = NULL;
	free ( set );
	set = NULL;
}

/*************************************************
Function:
    free_Sets2
Description:
    Free all the kmersets.
Input:
    1. sets:        the array of the kmerset
    2. num:     the number of kmerset
Output:
    None.
Return:
    None.
*************************************************/
PUBLIC_FUNC  void free_Sets2 ( KmerSet2 ** sets, int num )
{
	int i;

	for ( i = 0; i < num; ++i )
	{
		free_kmerset2 ( sets[i] );
		sets[i] = NULL;
	}

	free ( ( void * ) sets );
	sets = NULL;
}

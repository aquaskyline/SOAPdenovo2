/*
 * newhash.c
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

#ifdef MER127
static const kmer_t empty_kmer = { {0, 0, 0, 0}, 0, 0, 0, 0, 0, 1, 0, 0 };

static inline ubyte8 modular ( KmerSet * set, Kmer seq )
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
static const kmer_t empty_kmer = { {0, 0}, 0, 0, 0, 0, 0, 1, 0, 0 };
static inline ubyte8 modular ( KmerSet * set, Kmer seq )
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
    update_kmer
Description:
    Updates the coverage of the kmer.
Input:
    1. mer:     the current kmer
    2. left:        the upstream character
    3. right:       the downstream character
Output:
    None.
Return:
    None.
*************************************************/
static inline void update_kmer ( kmer_t * mer, ubyte left, ubyte right )
{
	ubyte4 cov;

	if ( left < 4 )
	{
		cov = get_kmer_left_cov ( *mer, left );

		if ( cov < MAX_KMER_COV )
		{
			set_kmer_left_cov ( *mer, left, cov + 1 );
		}
	}

	if ( right < 4 )
	{
		cov = get_kmer_right_cov ( *mer, right );

		if ( cov < MAX_KMER_COV )
		{
			set_kmer_right_cov ( *mer, right, cov + 1 );
		}
	}
}

/*************************************************
Function:
    set_new_kmer
Description:
    Initializes the new kmer.
Input:
    1. mer:     the point of the new kmer
    2. seq:     the sequence of the kmer
    3. left:        the upstream character
    4. right:       the downstream character
Output:
    None.
Return:
    None.
*************************************************/
static inline void set_new_kmer ( kmer_t * mer, Kmer seq, ubyte left, ubyte right )
{
	*mer = empty_kmer;
	set_kmer_seq ( *mer, seq );

	if ( left < 4 )
	{
		set_kmer_left_cov ( *mer, left, 1 );
	}

	if ( right < 4 )
	{
		set_kmer_right_cov ( *mer, right, 1 );
	}
}

static inline int is_prime_kh ( ubyte8 num )
{
	ubyte8 i, max;

	if ( num < 4 )
	{
		return 1;
	}

	if ( num % 2 == 0 )
	{
		return 0;
	}

	max = ( ubyte8 ) sqrt ( ( float ) num );

	for ( i = 3; i < max; i += 2 )
	{
		if ( num % i == 0 )
		{
			return 0;
		}
	}

	return 1;
}

static inline ubyte8 find_next_prime_kh ( ubyte8 num )
{
	if ( num % 2 == 0 )
	{
		num++;
	}

	while ( 1 )
	{
		if ( is_prime_kh ( num ) )
		{
			return num;
		}

		num += 2;
	}
}

/*************************************************
Function:
    init_kmerset
Description:
    Initializes the kmerset
Input:
    1. init_size:       the initial size of the kmerset
    2.load_factor:      load factor of hash
Output:
    None.
Return:
    The new kmer hash.
*************************************************/
PUBLIC_FUNC KmerSet * init_kmerset ( ubyte8 init_size, float load_factor )
{
	KmerSet * set;

	if ( init_size < 3 )
	{
		init_size = 3;
	}
	else
	{
		init_size = find_next_prime_kh ( init_size );
	}

	set = ( KmerSet * ) malloc ( sizeof ( KmerSet ) );
	set->size = init_size;
	set->count = 0;
	set->max = set->size * load_factor;

	if ( load_factor <= 0 )
	{
		load_factor = 0.25f;
	}
	else if ( load_factor >= 1 )
	{
		load_factor = 0.75f;
	}

	set->load_factor = load_factor;
	set->iter_ptr = 0;
	set->array = calloc ( set->size, sizeof ( kmer_t ) );
	set->flags = malloc ( ( set->size + 15 ) / 16 * 4 );
	memset ( set->flags, 0x55, ( set->size + 15 ) / 16 * 4 );
	return set;
}

PROTECTED_FUNC static inline ubyte8 get_kmerset ( KmerSet * set, Kmer seq )
{
	ubyte8 hc;
	//  U256b temp;
	//  temp = Kmer2int256(seq);
	//  hc = temp.low % set->size;
	//  hc = modular (set, seq);
	//  /*
#ifdef MER127
	hc = modular ( set, seq );
#else
	__uint128_t temp;
	temp = Kmer2int128 ( seq );
	hc = temp % set->size;
#endif
	//  */

	while ( 1 )
	{
		if ( is_kmer_entity_null ( set->flags, hc ) )
		{
			return hc;
		}
		else
		{
			if ( KmerEqual ( get_kmer_seq ( set->array[hc] ), seq ) )
			{
				return hc;
			}
		}

		hc++;

		if ( hc == set->size )
		{
			hc = 0;
		}
	}

	return set->size;
}

PUBLIC_FUNC int search_kmerset ( KmerSet * set, Kmer seq, kmer_t ** rs )
{
	ubyte8 hc;
	//  U256b temp;
	//  temp = Kmer2int256(seq);
	//  hc = temp.low % set->size;
	//  hc = modular (set, seq);
	//  /*
#ifdef MER127
	hc = modular ( set, seq );
#else
	__uint128_t temp;
	temp = Kmer2int128 ( seq );
	hc = temp % set->size;
#endif
	//  */

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

		hc++;

		if ( hc == set->size )
		{
			hc = 0;
		}
	}

	return 0;
}

PUBLIC_FUNC static inline int exists_kmerset ( KmerSet * set, Kmer seq )
{
	ubyte8 idx;
	idx = get_kmerset ( set, seq );
	return !is_kmer_entity_null ( set->flags, idx );
}

/*************************************************
Function:
    encap_kmerset
Description:
    Enlarges the kmerset if necessary.
Input:
    1. set :        the hash of kmerset
    2. num:     the element number add to kmerset
Output:
    None.
Return:
    None.
*************************************************/
PROTECTED_FUNC static inline void encap_kmerset ( KmerSet * set, ubyte8 num )
{
	ubyte4 * flags, *f;
	ubyte8 i, n, size, hc;
	kmer_t key, tmp;

	if ( set->count + num <= set->max )
	{
		return;
	}

	n = set->size;

	if ( initKmerSetSize != 0 )
	{
		if ( set->load_factor < 0.88 )
		{
			set->load_factor = 0.88;
			set->max = set->size * set->load_factor;
			return;
		}
		else
		{
			fprintf ( stderr, "-- Static memory pool exploded, please define a larger value. --\n" );
			abort();
		}
	}

	do
	{
		if ( n < 0xFFFFFFFU )
		{
			n <<= 1;
		}
		else
		{
			n += 0xFFFFFFU;
		}

		n = find_next_prime_kh ( n );
	}
	while ( n * set->load_factor < set->count + num );

	set->array = realloc ( set->array, n * sizeof ( kmer_t ) );
	//printf("Allocate Mem %lld(%d*%lld*%d)bytes\n",thrd_num*n*sizeof(kmer_t),thrd_num,n,sizeof(kmer_t));

	if ( set->array == NULL )
	{
		fprintf ( stderr, "-- Out of memory --\n" );
		abort ();
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

	//  U256b temp;
	for ( i = 0; i < size; i++ )
	{
		if ( !exists_kmer_entity ( flags, i ) )
		{
			continue;
		}

		key = set->array[i];
		set_kmer_entity_del ( flags, i );

		while ( 1 )
		{
			//          temp = Kmer2int256(get_kmer_seq(key));
			//          hc = temp.low % set->size;
			//          hc = modular (set, get_kmer_seq (key));
			//          /*
#ifdef MER127
			hc = modular ( set, get_kmer_seq ( key ) );
#else
			temp = Kmer2int128 ( get_kmer_seq ( key ) );
			hc = temp % set->size;
#endif
			//          */

			while ( !is_kmer_entity_null ( set->flags, hc ) )
			{
				hc++;

				if ( hc == set->size )
				{
					hc = 0;
				}
			}

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
    put_kmerset
Description:
    Put a kmer into the kmerset.
Input:
    1. set :        the hash of kmerset
    2. seq:     the sequence of the kmer
    3. left:        the upstream charactor of the kmer
    4. right:       the downstream charactor of the kmer
    5. kmer_p:  the pointer to the kmer
Output:
    None.
Return:
    0 if it failed to put kmer into kmerset.
*************************************************/
PUBLIC_FUNC int put_kmerset ( KmerSet * set, Kmer seq, ubyte left, ubyte right, kmer_t ** kmer_p )
{
	ubyte8 hc;

	if ( set->count + 1 > set->max )
	{
		encap_kmerset ( set, 1 );
	}

	//  U256b temp;
	//  temp = Kmer2int256(seq);
	//  hc = temp.low % set->size;
	//  hc = modular (set, seq);
	//  /*
#ifdef MER127
	hc = modular ( set, seq );
#else
	__uint128_t temp;
	temp = Kmer2int128 ( seq );
	hc = temp % set->size;
#endif
	//  */

	do
	{
		if ( is_kmer_entity_null ( set->flags, hc ) )
		{
			clear_kmer_entity_null ( set->flags, hc );
			set_new_kmer ( set->array + hc, seq, left, right );
			set->count++;
			*kmer_p = set->array + hc;
			return 0;
		}
		else
		{
			if ( KmerEqual ( get_kmer_seq ( set->array[hc] ), seq ) )
			{
				update_kmer ( set->array + hc, left, right );
				set->array[hc].single = 0;
				*kmer_p = set->array + hc;
				return 1;
			}
		}

		hc++;

		if ( hc == set->size )
		{
			hc = 0;
		}
	}
	while ( 1 );

	*kmer_p = NULL;
	return 0;
}

/*************************************************
Function:
    count_kmerset
Description:
    1. Returns the kmer number of the kmerset.
Input:
    1. set:     the kmer hash
Output:
    None.
Return:
    The kmer number of the kmerset.
*************************************************/
PUBLIC_FUNC byte8 count_kmerset ( KmerSet * set )
{
	return set->count;
}

/*************************************************
Function:
    reset_iter_kmerset
Description:
    Resets the iterator of kmerset.
Input:
    1. set:     the kmer hash
Output:
    None.
Return:
    None.
*************************************************/
PUBLIC_FUNC static inline void reset_iter_kmerset ( KmerSet * set )
{
	set->iter_ptr = 0;
}

PUBLIC_FUNC static inline ubyte8 iter_kmerset ( KmerSet * set, kmer_t ** rs )
{
	while ( set->iter_ptr < set->size )
	{
		if ( !is_kmer_entity_null ( set->flags, set->iter_ptr ) )
		{
			*rs = set->array + set->iter_ptr;
			set->iter_ptr++;
			return 1;
		}

		set->iter_ptr++;
	}

	return 0;
}

PUBLIC_FUNC void free_kmerset ( KmerSet * set )
{
	free ( set->array );
	free ( set->flags );
	free ( set );
}

/*************************************************
Function:
    free_Sets
Description:
    Releases all the kmerset.
Input:
    1. set:     the kmer hash
    2. num:     the number of kmersets
Output:
    None.
Return:
    None.
*************************************************/
PUBLIC_FUNC void free_Sets ( KmerSet ** sets, int num )
{
	int i;

	for ( i = 0; i < num; i++ )
	{
		free_kmerset ( sets[i] );
	}

	free ( ( void * ) sets );
}

/*************************************************
Function:
    count_branch2prev
Description:
    Calculates the number of the upstream branches of a kmer.
Input:
    1. node:        the kmer
Output:
    None.
Return:
    The number of the upstream branches.
*************************************************/
int count_branch2prev ( kmer_t * node )
{
	int num = 0, i;

	for ( i = 0; i < 4; i++ )
	{
		if ( get_kmer_left_cov ( *node, i ) > 0 )
		{
			num++;
		}
	}

	return num;
}

/*************************************************
Function:
    count_branch2next
Description:
    Calculates the number of the downstream branches of a kmer.
Input:
    1. node:        the kmer
Output:
    None.
Return:
    The number of the downstream branches.
*************************************************/
int count_branch2next ( kmer_t * node )
{
	int num = 0, i;

	for ( i = 0; i < 4; i++ )
	{
		if ( get_kmer_right_cov ( *node, i ) > 0 )
		{
			num++;
		}
	}

	return num;
}

/*************************************************
Function:
    dislink2prevUncertain
Description:
    Removes the relation between a kmer and an upstream base.
Input:
    1. node:        the current kmer
    2. ch:      the upstream base
    3. smaller: the direction of the kmer
Output:
    None.
Return:
    None.
*************************************************/
void dislink2prevUncertain ( kmer_t * node, char ch, boolean smaller )
{
	if ( smaller )
	{
		set_kmer_left_cov ( *node, ch, 0 );
	}
	else
	{
		set_kmer_right_cov ( *node, int_comp ( ch ), 0 );
	}
}

/*************************************************
Function:
    dislink2nextUncertain
Description:
    Removes the relation between a kmer and a downstream base.
Input:
    1. node:        the current kmer
    2. ch:      the downstream base
    3. smaller: the direction of the kmer
Output:
    None.
Return:
    None.
*************************************************/
void dislink2nextUncertain ( kmer_t * node, char ch, boolean smaller )
{
	if ( smaller )
	{
		set_kmer_right_cov ( *node, ch, 0 );
	}
	else
	{
		set_kmer_left_cov ( *node, int_comp ( ch ), 0 );
	}
}

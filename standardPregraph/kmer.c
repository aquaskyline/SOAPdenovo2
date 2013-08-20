/*
 * kmer.c
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

#ifdef MER127
void PrintKmer ( FILE * fp, Kmer kmer )
{
	fprintf ( fp, "%llx %llx %llx %llx", kmer.high1, kmer.low1, kmer.high2, kmer.low2 );
}

/*************************************************
Function:
    KmerSmaller
Description:
    Checks if kmer1 is smaller than kmer2. if no, returns 0.
Input:
    1. kmer1:       the first kmer
    2. kmer2:       the second kmer
Output:
    None.
Return:
    1 if kmer1 is smaller than kmer2.
*************************************************/
boolean KmerSmaller ( Kmer kmer1, Kmer kmer2 )
{
	if ( kmer1.high1 != kmer2.high1 )
	{
		return ( kmer1.high1 < kmer2.high1 );
	}
	else
	{
		if ( kmer1.low1 != kmer2.low1 )
		{
			return ( kmer1.low1 < kmer2.low1 );
		}
		else
		{
			if ( kmer1.high2 != kmer2.high2 )
			{
				return ( kmer1.high2 < kmer2.high2 );
			}
			else
			{
				return ( kmer1.low2 < kmer2.low2 );
			}
		}
	}
}

/*************************************************
Function:
    KmerLarger
Description:
    Checks if kmer1 is larger than kmer2. if no, returns 0.
Input:
    1. kmer1:       the first kmer
    2. kmer2:       the second kmer
Output:
    None.
Return:
    1 if kmer1 is larger than kmer2.
*************************************************/
boolean KmerLarger ( Kmer kmer1, Kmer kmer2 )
{
	if ( kmer1.high1 != kmer2.high1 )
	{
		return ( kmer1.high1 > kmer2.high1 );
	}
	else
	{
		if ( kmer1.low1 != kmer2.low1 )
		{
			return ( kmer1.low1 > kmer2.low1 );
		}
		else
		{
			if ( kmer1.high2 != kmer2.high2 )
			{
				return ( kmer1.high2 > kmer2.high2 );
			}
			else
			{
				return ( kmer1.low2 > kmer2.low2 );
			}
		}
	}
}

/*************************************************
Function:
    KmerEqual
Description:
    Checks if kmer1 is same as kmer2. if no, returns 0.
Input:
    1. kmer1:       the first kmer
    2. kmer2:       the second kmer
Output:
    None.
Return:
    1 if if kmer1 is same as kmer2.
*************************************************/
boolean KmerEqual ( Kmer kmer1, Kmer kmer2 )
{
	if ( kmer1.low2 != kmer2.low2 || kmer1.high2 != kmer2.high2 || kmer1.low1 != kmer2.low1 || kmer1.high1 != kmer2.high1 )
	{
		return 0;
	}
	else
	{
		return 1;
	}
}

/*************************************************
Function:
    KmerAnd
Description:
    kmer1 & kmer2
Input:
    1. kmer1:       the first kmer
    2. kmer2:       the second kmer
Output:
    None.
Return:
    kmer1 & kmer2.
*************************************************/
Kmer KmerAnd ( Kmer kmer1, Kmer kmer2 )
{
	kmer1.high1 &= kmer2.high1;
	kmer1.low1 &= kmer2.low1;
	kmer1.high2 &= kmer2.high2;
	kmer1.low2 &= kmer2.low2;
	return kmer1;
}

/*************************************************
Function:
    KmerLeftBitMoveBy2
Description:
    word <<= 2.
Input:
    1. word:        the kmer
Output:
    None.
Return:
    word <<= 2.
*************************************************/
Kmer KmerLeftBitMoveBy2 ( Kmer word )
{
	word.high1 = ( word.high1 << 2 ) | ( word.low1 >> 62 );
	word.low1 = ( word.low1 << 2 ) | ( word.high2 >> 62 );
	word.high2 = ( word.high2 << 2 ) | ( word.low2 >> 62 );
	word.low2 <<= 2;
	return word;
}

/*************************************************
Function:
    KmerRightBitMoveBy2
Description:
    word >>= 2.
Input:
    1. word:        the kmer
Output:
    None.
Return:
    word >>= 2.
*************************************************/
Kmer KmerRightBitMoveBy2 ( Kmer word )
{
	word.low2 = ( word.low2 >> 2 ) | ( word.high2 & 0x3 ) << 62;
	word.high2 = ( word.high2 >> 2 ) | ( word.low1 & 0x3 ) << 62;
	word.low1 = ( word.low1 >> 2 ) | ( word.high1 & 0x3 ) << 62;
	word.high1 >>= 2;
	return word;
}

/*************************************************
Function:
    KmerPlus
Description:
    Appends a base to kmer to get (k+1)mer.
Input:
    1. prev:        the kmer
    2. ch:      the downstream
Output:
    None.
Return:
    (k+1)mer.
*************************************************/
Kmer KmerPlus ( Kmer prev, char ch )
{
	Kmer word = KmerLeftBitMoveBy2 ( prev );
	word.low2 |= ch;
	return word;
}

/*************************************************
Function:
    nextKmer
Description:
    Combines the last (k-1) bases of a kmer with another base to get new kmer.
Input:
    1. prev:        the kmer
    2. ch:      the appened base
Output:
    None.
Return:
    The new kmer.
*************************************************/
Kmer nextKmer ( Kmer prev, char ch )
{
	Kmer word = KmerLeftBitMoveBy2 ( prev );
	word = KmerAnd ( word, WORDFILTER );
	word.low2 |= ch;
	return word;
}

/*************************************************
Function:
    prevKmer
Description:
    Combines a base with the first (k-1) bases of a kmer to get new kmer.
Input:
    1. next:        the kmer
    2. ch:      the first base of new kmer
Output:
    None.
Return:
    The new kmer.
*************************************************/
Kmer prevKmer ( Kmer next, char ch )
{
	Kmer word = KmerRightBitMoveBy2 ( next );

	switch ( overlaplen )
	{
		case 1 ... 32:
			word.low2 |= ( ( ( ubyte8 ) ch ) << 2 * ( overlaplen - 1 ) );
			break;
		case 33 ... 64:
			word.high2 |= ( ( ubyte8 ) ch ) << ( 2 * ( overlaplen - 1 ) - 64 );
			break;
		case 65 ... 96:
			word.low1 |= ( ( ubyte8 ) ch ) << ( 2 * ( overlaplen - 1 ) - 128 );
			break;
		case 97 ... 128:
			word.high1 |= ( ( ubyte8 ) ch ) << ( 2 * ( overlaplen - 1 ) - 192 );
			break;
	}

	return word;
}

/*************************************************
Function:
    lastCharInKmer
Description:
    Gets the last char of a kmer, and returns it.
Input:
    1. kmer:        the kmer
Output:
    None.
Return:
    The last char of kmer.
*************************************************/
char lastCharInKmer ( Kmer kmer )
{
	return ( char ) ( kmer.low2 & 0x3 );
}

/*************************************************
Function:
    firstCharInKmer
Description:
    Gets the first char of a kmer, and returns it.
Input:
    1. kmer:        the kmer
Output:
    None.
Return:
    The first char of kmer.
*************************************************/
char firstCharInKmer ( Kmer kmer )
{
	switch ( overlaplen )
	{
		case 1 ... 32:
			kmer.low2 >>= 2 * ( overlaplen - 1 );
			return kmer.low2;   // & 3;
		case 33 ... 64:
			kmer.high2 >>= 2 * ( overlaplen - 1 ) - 64;
			return kmer.high2;  // & 3;
		case 65 ... 96:
			kmer.low1 >>= 2 * ( overlaplen - 1 ) - 128;
			return kmer.low1;
		case 97 ... 128:
			kmer.high1 >>= 2 * ( overlaplen - 1 ) - 192;
			return kmer.high1;
	}
}

/*************************************************
Function:
    createFilter
Description:
    Creates filter for masking.
Input:
    1. overlaplen:      the kmer value
Output:
    None.
Return:
    The filter kmer.
*************************************************/
Kmer createFilter ( int overlaplen )
{
	Kmer word;
	word.high1 = word.low1 = word.high2 = word.low2 = 0;

	switch ( overlaplen )
	{
		case 1 ... 31:
			word.low2 = ( ( ( ubyte8 ) 1 ) << ( 2 * overlaplen ) ) - 1;
			break;
		case 32 ... 63:
			word.low2 = ~word.low2;
			word.high2 = ( ( ( ubyte8 ) 1 ) << ( 2 * overlaplen - 64 ) ) - 1;
			break;
		case 64 ... 95:
			word.high2 = word.low2 = ~word.low2;
			word.low1 = ( ( ( ubyte8 ) 1 ) << ( 2 * overlaplen - 128 ) ) - 1;
			break;
		case 96 ... 127:
			word.low1 = word.high2 = word.low2 = ~word.low2;
			word.high1 = ( ( ( ubyte8 ) 1 ) << ( 2 * overlaplen - 192 ) ) - 1;
			break;
	}

	return word;
}

/*************************************************
Function:
    KmerRightBitMove
Description:
    Makes right shifts to a kmer by specified bits and returns the new kmer.
Input:
    1. word:        the kmer
    2. dis:     the shifted bits
Output:
    None.
Return:
    The new kmer.
*************************************************/
Kmer KmerRightBitMove ( Kmer word, int dis )
{
	ubyte8 mask;

	switch ( dis )
	{
		case 0 ... 63:
			mask = ( ( ( ubyte8 ) 1 ) << dis ) - 1;
			word.low2 = ( word.low2 >> dis ) | ( word.high2 & mask ) << ( 64 - dis );
			word.high2 = ( word.high2 >> dis ) | ( word.low1 & mask ) << ( 64 - dis );
			word.low1 = ( word.low1 >> dis ) | ( word.high1 & mask ) << ( 64 - dis );
			word.high1 >>= dis;
			return word;
		case 64 ... 127:
			mask = ( ( ( ubyte8 ) 1 ) << ( dis - 64 ) ) - 1;
			word.low2 = word.high2 >> ( dis - 64 ) | ( word.low1 & mask ) << ( 128 - dis );
			word.high2 = word.low1 >> ( dis - 64 ) | ( word.high1 & mask ) << ( 128 - dis );
			word.low1 = word.high1 >> ( dis - 64 );
			word.high1 = 0;
			return word;
		case 128 ... 191:
			mask = ( ( ( ubyte8 ) 1 ) << ( dis - 128 ) ) - 1;
			word.low2 = word.low1 >> ( dis - 128 ) | ( word.high1 & mask ) << ( 192 - dis );
			word.high2 = word.high1 >> ( dis - 128 );
			word.high1 = word.low1 = 0;
			return word;
		case 192 ... 255:
			word.low2 = word.high1 >> ( dis - 192 );
			word.high1 = word.low1 = word.high2 = 0;
			return word;
	}
}

/*************************************************
Function:
    printKmerSeq
Description:
    Prints seq of kmer.
Input:
    1. fp:      output file
    2. kmer:        the output kmer
Output:
    None.
Return:
    None.
*************************************************/
void printKmerSeq ( FILE * fp, Kmer kmer )
{
	int i, bit1, bit2, bit3, bit4;
	bit4 = bit3 = bit2 = bit1 = 0;
	char kmerSeq[128];

	switch ( overlaplen )
	{
		case 1 ... 31:
			bit4 = overlaplen;
			break;
		case 32 ... 63:
			bit4 = 32;
			bit3 = overlaplen - 32;
			break;
		case 64 ... 95:
			bit4 = bit3 = 32;
			bit2 = overlaplen - 64;
			break;
		case 96 ... 127:
			bit4 = bit3 = bit2 = 32;
			bit1 = overlaplen - 96;
			break;
	}

	for ( i = bit1 - 1; i >= 0; i-- )
	{
		kmerSeq[i] = kmer.high1 & 0x3;
		kmer.high1 >>= 2;
	}

	for ( i = bit2 - 1; i >= 0; i-- )
	{
		kmerSeq[i + bit1] = kmer.low1 & 0x3;
		kmer.low1 >>= 2;
	}

	for ( i = bit3 - 1; i >= 0; i-- )
	{
		kmerSeq[i + bit1 + bit2] = kmer.high2 & 0x3;
		kmer.high2 >>= 2;
	}

	for ( i = bit4 - 1; i >= 0; i-- )
	{
		kmerSeq[i + bit1 + bit2 + bit3] = kmer.low2 & 0x3;
		kmer.low2 >>= 2;
	}

	for ( i = 0; i < overlaplen; i++ )
	{
		fprintf ( fp, "%c", int2base ( ( int ) kmerSeq[i] ) );
	}
}

void print_kmer ( FILE * fp, Kmer kmer, char c )
{
	fprintf ( fp, "%llx %llx %llx %llx", kmer.high1, kmer.low1, kmer.high2, kmer.low2 );
	fprintf ( fp, "%c", c );
}

void print_kmer_gz ( gzFile * fp, Kmer kmer, char c )
{
	gzprintf ( fp, "%llx %llx %llx %llx", kmer.high1, kmer.low1, kmer.high2, kmer.low2 );
	gzprintf ( fp, "%c", c );
}

static const ubyte2 BitReverseTable[65536] =
{
# define R2(n)      n,     n + 1*16384,    n + 2*16384,    n + 3*16384
# define R4(n)   R2(n), R2(n + 1*4096), R2(n + 2*4096), R2(n + 3*4096)
# define R6(n)   R4(n), R4(n + 1*1024), R4(n + 2*1024), R4(n + 3*1024)
# define R8(n)   R6(n), R6(n + 1*256 ), R6(n + 2*256 ), R6(n + 3*256 )
# define R10(n)  R8(n), R8(n + 1*64  ), R8(n + 2*64  ), R8(n + 3*64  )
# define R12(n)  R10(n),R10(n + 1*16), R10(n + 2*16 ), R10(n + 3*16  )
# define R14(n)  R12(n),R12(n + 1*4 ), R12(n + 2*4  ), R12(n + 3*4  )
	R14 ( 0 ), R14 ( 1 ), R14 ( 2 ), R14 ( 3 )
};

/*************************************************
Function:
    fastReverseComp
Description:
    Gets the reverse complement of a kmer.
Input:
    1. seq:         the kmer
    2. seq_size:        the length of the kmer
Output:
    None.
Return:
    The reverse complement of kmer.
*************************************************/
static Kmer fastReverseComp ( Kmer seq, char seq_size )
{
	seq.low2 ^= 0xAAAAAAAAAAAAAAAALLU;
	seq.low2 = ( ( ubyte8 ) BitReverseTable[seq.low2 & 0xffff] << 48 ) |
	           ( ( ubyte8 ) BitReverseTable[ ( seq.low2 >> 16 ) & 0xffff] << 32 ) | ( ( ubyte8 ) BitReverseTable[ ( seq.low2 >> 32 ) & 0xffff] << 16 ) | ( ( ubyte8 ) BitReverseTable[ ( seq.low2 >> 48 ) & 0xffff] );

	if ( seq_size < 32 )
	{
		seq.low2 >>= ( 64 - ( seq_size << 1 ) );
		return seq;
	}

	seq.high2 ^= 0xAAAAAAAAAAAAAAAALLU;
	seq.high2 = ( ( ubyte8 ) BitReverseTable[seq.high2 & 0xffff] << 48 ) |
	            ( ( ubyte8 ) BitReverseTable[ ( seq.high2 >> 16 ) & 0xffff] << 32 ) | ( ( ubyte8 ) BitReverseTable[ ( seq.high2 >> 32 ) & 0xffff] << 16 ) | ( ( ubyte8 ) BitReverseTable[ ( seq.high2 >> 48 ) & 0xffff] );

	if ( seq_size < 64 )
	{
		seq.high2 = seq.high2 ^ seq.low2;
		seq.low2 = seq.high2 ^ seq.low2;
		seq.high2 = seq.high2 ^ seq.low2;
		seq = KmerRightBitMove ( seq, 128 - ( seq_size << 1 ) );
		return seq;
	}

	seq.low1 ^= 0xAAAAAAAAAAAAAAAALLU;
	seq.low1 = ( ( ubyte8 ) BitReverseTable[seq.low1 & 0xffff] << 48 ) |
	           ( ( ubyte8 ) BitReverseTable[ ( seq.low1 >> 16 ) & 0xffff] << 32 ) | ( ( ubyte8 ) BitReverseTable[ ( seq.low1 >> 32 ) & 0xffff] << 16 ) | ( ( ubyte8 ) BitReverseTable[ ( seq.low1 >> 48 ) & 0xffff] );

	if ( seq_size < 96 )
	{
		seq.low1 = seq.low1 ^ seq.low2;
		seq.low2 = seq.low1 ^ seq.low2;
		seq.low1 = seq.low1 ^ seq.low2;
		seq = KmerRightBitMove ( seq, 192 - ( seq_size << 1 ) );
		return seq;
	}

	seq.high1 ^= 0xAAAAAAAAAAAAAAAALLU;
	seq.high1 = ( ( ubyte8 ) BitReverseTable[seq.high1 & 0xffff] << 48 ) |
	            ( ( ubyte8 ) BitReverseTable[ ( seq.high1 >> 16 ) & 0xffff] << 32 ) | ( ( ubyte8 ) BitReverseTable[ ( seq.high1 >> 32 ) & 0xffff] << 16 ) | ( ( ubyte8 ) BitReverseTable[ ( seq.high1 >> 48 ) & 0xffff] );
	seq.low1 = seq.low1 ^ seq.high2;
	seq.high2 = seq.low1 ^ seq.high2;
	seq.low1 = seq.low1 ^ seq.high2;
	seq.low2 = seq.low2 ^ seq.high1;
	seq.high1 = seq.low2 ^ seq.high1;
	seq.low2 = seq.low2 ^ seq.high1;
	seq = KmerRightBitMove ( seq, 256 - ( seq_size << 1 ) );
	return seq;
}

Kmer reverseComplementVerbose ( Kmer word, int overlap )
{
	return fastReverseComp ( word, overlap );
}

Kmer reverseComplement ( Kmer word, int overlap )
{
	return fastReverseComp ( word, overlap );
}

#else
void PrintKmer ( FILE * fp, Kmer kmer )
{
	fprintf ( fp, "%llx %llx", kmer.high, kmer.low );
}


__uint128_t Kmer2int128 ( Kmer seq )
{
	__uint128_t temp;
	temp = seq.high;
	temp <<= 64;
	temp |= seq.low;
	return temp;
}
boolean KmerSmaller ( Kmer kmer1, Kmer kmer2 )
{
	if ( kmer1.high < kmer2.high )
	{
		return 1;
	}
	else if ( kmer1.high == kmer2.high )
	{
		if ( kmer1.low < kmer2.low )
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
	else
	{
		return 0;
	}
}

boolean KmerLarger ( Kmer kmer1, Kmer kmer2 )
{
	if ( kmer1.high > kmer2.high )
	{
		return 1;
	}
	else if ( kmer1.high == kmer2.high )
	{
		if ( kmer1.low > kmer2.low )
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
	else
	{
		return 0;
	}
}
boolean KmerEqual ( Kmer kmer1, Kmer kmer2 )
{
	if ( kmer1.high == kmer2.high && kmer1.low == kmer2.low )
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

Kmer KmerAnd ( Kmer kmer1, Kmer kmer2 )
{
	kmer1.high &= kmer2.high;
	kmer1.low &= kmer2.low;
	return kmer1;
}

Kmer KmerLeftBitMoveBy2 ( Kmer word )
{
	ubyte8 temp = word.low >> 62;
	word.high <<= 2;
	word.high |= temp;
	word.low <<= 2;
	return word;
}

Kmer KmerRightBitMoveBy2 ( Kmer word )
{
	ubyte8 temp = ( word.high & 0x3 ) << 62;
	word.high >>= 2;
	word.low >>= 2;
	word.low |= temp;
	return word;
}

Kmer KmerPlus ( Kmer prev, char ch )
{
	Kmer word = KmerLeftBitMoveBy2 ( prev );
	word.low |= ch;
	return word;
}
Kmer nextKmer ( Kmer prev, char ch )
{
	Kmer word = KmerLeftBitMoveBy2 ( prev );
	word = KmerAnd ( word, WORDFILTER );
	word.low |= ch;
	return word;
}

Kmer prevKmer ( Kmer next, char ch )
{
	Kmer word = KmerRightBitMoveBy2 ( next );

	if ( 2 * ( overlaplen - 1 ) < 64 )
	{
		word.low |= ( ( ( ubyte8 ) ch ) << 2 * ( overlaplen - 1 ) );
	}
	else
	{
		word.high |= ( ( ubyte8 ) ch ) << ( 2 * ( overlaplen - 1 ) - 64 );
	}

	return word;
}

char lastCharInKmer ( Kmer kmer )
{
	return ( char ) ( kmer.low & 0x3 );
}
char firstCharInKmer ( Kmer kmer )
{
	if ( 2 * ( overlaplen - 1 ) < 64 )
	{
		kmer.low >>= 2 * ( overlaplen - 1 );
		return kmer.low;    // & 3;
	}
	else
	{
		kmer.high >>= 2 * ( overlaplen - 1 ) - 64;
		return kmer.high;   // & 3;
	}
}

Kmer createFilter ( int overlaplen )
{
	Kmer word;
	word.high = word.low = 0;

	if ( 2 * overlaplen < 64 )
	{
		word.low = ( ( ( ubyte8 ) 1 ) << ( 2 * overlaplen ) ) - 1;
	}
	else
	{
		word.low = ~word.low;

		if ( 2 * overlaplen > 64 )
		{
			word.high = ( ( ( ubyte8 ) 1 ) << ( 2 * overlaplen - 64 ) ) - 1;
		}
	}

	return word;
}

Kmer KmerRightBitMove ( Kmer word, int dis )
{
	if ( dis < 64 )
	{
		ubyte8 mask = ( ( ( ubyte8 ) 1 ) << dis ) - 1;
		ubyte8 temp = ( word.high & mask ) << ( 64 - dis );
		word.high >>= dis;
		word.low >>= dis;
		word.low |= temp;
		return word;
	}

	word.high >>= ( dis - 64 );
	word.low = word.high;
	word.high = 0;
	return word;
}


void printKmerSeq ( FILE * fp, Kmer kmer )
{
	int i, bit1, bit2;
	char ch;
	char kmerSeq[64];
	bit2 = overlaplen > 32 ? 32 : overlaplen;
	bit1 = overlaplen > 32 ? overlaplen - 32 : 0;

	for ( i = bit1 - 1; i >= 0; i-- )
	{
		ch = kmer.high & 0x3;
		kmer.high >>= 2;
		kmerSeq[i] = ch;
	}

	for ( i = bit2 - 1; i >= 0; i-- )
	{
		ch = kmer.low & 0x3;
		kmer.low >>= 2;
		kmerSeq[i + bit1] = ch;
	}

	for ( i = 0; i < overlaplen; i++ )
	{
		fprintf ( fp, "%c", int2base ( ( int ) kmerSeq[i] ) );
	}
}

void print_kmer ( FILE * fp, Kmer kmer, char c )
{
	fprintf ( fp, "%llx %llx", kmer.high, kmer.low );
	fprintf ( fp, "%c", c );
}

void print_kmer_gz ( gzFile * fp, Kmer kmer, char c )
{
	gzprintf ( fp, "%llx %llx", kmer.high, kmer.low );
	gzprintf ( fp, "%c", c );
}

static Kmer fastReverseComp ( Kmer seq, char seq_size )
{
	seq.low ^= 0xAAAAAAAAAAAAAAAALLU;
	seq.low = ( ( seq.low & 0x3333333333333333LLU ) << 2 ) | ( ( seq.low & 0xCCCCCCCCCCCCCCCCLLU ) >> 2 );
	seq.low = ( ( seq.low & 0x0F0F0F0F0F0F0F0FLLU ) << 4 ) | ( ( seq.low & 0xF0F0F0F0F0F0F0F0LLU ) >> 4 );
	seq.low = ( ( seq.low & 0x00FF00FF00FF00FFLLU ) << 8 ) | ( ( seq.low & 0xFF00FF00FF00FF00LLU ) >> 8 );
	seq.low = ( ( seq.low & 0x0000FFFF0000FFFFLLU ) << 16 ) | ( ( seq.low & 0xFFFF0000FFFF0000LLU ) >> 16 );
	seq.low = ( ( seq.low & 0x00000000FFFFFFFFLLU ) << 32 ) | ( ( seq.low & 0xFFFFFFFF00000000LLU ) >> 32 );

	if ( seq_size < 32 )
	{
		seq.low >>= ( 64 - ( seq_size << 1 ) );
		return seq;
	}

	seq.high ^= 0xAAAAAAAAAAAAAAAALLU;
	seq.high = ( ( seq.high & 0x3333333333333333LLU ) << 2 ) | ( ( seq.high & 0xCCCCCCCCCCCCCCCCLLU ) >> 2 );
	seq.high = ( ( seq.high & 0x0F0F0F0F0F0F0F0FLLU ) << 4 ) | ( ( seq.high & 0xF0F0F0F0F0F0F0F0LLU ) >> 4 );
	seq.high = ( ( seq.high & 0x00FF00FF00FF00FFLLU ) << 8 ) | ( ( seq.high & 0xFF00FF00FF00FF00LLU ) >> 8 );
	seq.high = ( ( seq.high & 0x0000FFFF0000FFFFLLU ) << 16 ) | ( ( seq.high & 0xFFFF0000FFFF0000LLU ) >> 16 );
	seq.high = ( ( seq.high & 0x00000000FFFFFFFFLLU ) << 32 ) | ( ( seq.high & 0xFFFFFFFF00000000LLU ) >> 32 );
	ubyte8 temp = seq.high;
	seq.high = seq.low;
	seq.low = temp;
	seq = KmerRightBitMove ( seq, 128 - ( seq_size << 1 ) );
	return seq;
}

Kmer reverseComplementVerbose ( Kmer word, int overlap )
{
	return fastReverseComp ( word, overlap );
}

Kmer reverseComplement ( Kmer word, int overlap )
{
	return fastReverseComp ( word, overlap );
}

#endif

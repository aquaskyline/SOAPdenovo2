/*
 * inc/seq_util.h
 *
 * This file is part of SOAPdenovo.
 *
 * Part of this file is refered and modified from SparseAssembler
 * (See <http://sourceforge.net/projects/sparseassembler/>).
 *
 */

#ifndef _SEQ_UTIL_H
#define _SEQ_UTIL_H

#include "stdinc.h"
#include "core.h"


extern inline void Init_Read ( string &seq, struct read_t &read );

extern inline uint64_t *str2bitsarr ( const char *c_str, int len, uint64_t *b_str, int arr_sz );

extern inline char *bitsarr2str ( uint64_t *b_seq, int len, char *c_str, int arr_sz );

extern inline void get_sub_arr ( uint64_t *bitsarr_in, int bitsarr_len, int begin_pos, int sub_sz, uint64_t *bitsarr_out );

extern inline void L_shift_NC ( uint64_t *bitsarr, int shift_sz, int arr_sz );

extern inline void R_shift_NC ( uint64_t *bitsarr, int shift_sz, int arr_sz );

extern inline int uint64_t_cmp ( uint64_t *A, uint64_t *B, int Kmer_arr_sz );

extern inline uint64_t *get_rev_comp_seq_arr ( uint64_t *seq_arr, int seq_size, int arr_sz );

extern inline uint64_t get_rev_comp_seq ( uint64_t seq, int seq_size );

extern inline uint64_t MurmurHash64A ( const void *key, int len, unsigned int seed );




inline void Init_Read ( string &seq, struct read_t &read )
{
  read.readLen = ( int ) seq.size();
  int Read_arr_sz = read.readLen / 32 + 1;
  int rem = read.readLen % 32;

  if ( rem == 0 )
    {
      Read_arr_sz--;
    }

  str2bitsarr ( seq.c_str(), ( int ) seq.size(), read.read_bits, Read_arr_sz );
}


inline uint64_t *str2bitsarr ( const char *c_str, int len, uint64_t *b_str, int arr_sz )
{
  for ( int  k = 0; k < arr_sz; ++k )
    {
      b_str[k] = 0;
    }

  int arr_sz_needed = len / 32 + 1;
  int rem = len % 32;

  if ( rem == 0 )
    {
      arr_sz_needed--;
    }

  int beg_arr_idx = arr_sz - arr_sz_needed;

  if ( rem == 0 && arr_sz_needed > 0 )
    {
      rem = 32;
    }

  for ( int k = 0; k < len; k++ )
    {
      if ( rem == 0 )
        {
          beg_arr_idx++;
          rem = 32;
        }

      switch ( c_str[k] )
        {
        case ( 'A' ) :
        case ( 'a' ) :
        case ( '0' ) :
          b_str[beg_arr_idx] <<= 2;
          rem--;
          break;

        case ( 'C' ) :
        case ( 'c' ) :
        case ( '1' ) :
          b_str[beg_arr_idx] <<= 2;
          ++b_str[beg_arr_idx];
          rem--;
          break;

        case 'G':
        case 'g':
        case '2':
          b_str[beg_arr_idx] <<= 1;
          ++b_str[beg_arr_idx];
          b_str[beg_arr_idx] <<= 1;
          rem--;
          break;

        case 'T':
        case 't':
        case '3':
          b_str[beg_arr_idx] <<= 1;
          ++b_str[beg_arr_idx];
          b_str[beg_arr_idx] <<= 1;
          ++b_str[beg_arr_idx];
          rem--;
          break;

        default:
          return b_str;
        }
    }

  return b_str;
}

inline char *bitsarr2str ( uint64_t *b_seq, int len, char *c_str, int arr_sz )
{
  int tot_bits = arr_sz * 64;

  for ( int i = 0; i < len; ++i )
    {
      uint64_t temp, temp2[100];

      for ( int k = 0; k < arr_sz; ++k )
        {
          temp2[k] = b_seq[k];
        }

      L_shift_NC ( temp2, tot_bits - ( len - i ) * 2, arr_sz );
      R_shift_NC ( temp2, tot_bits - 2, arr_sz );
      temp = temp2[arr_sz - 1];

      switch ( temp )
        {
        case 0:
          c_str[i] = 'A';
          break;

        case 1:
          c_str[i] = 'C';
          break;

        case 2:
          c_str[i] = 'G';
          break;

        case 3:
          c_str[i] = 'T';
          break;
        }
    }

  c_str[len] = '\0';
  return c_str;
}



inline void get_sub_arr ( uint64_t *bitsarr_in, int bitsarr_len, int begin_pos, int sub_sz, uint64_t *bitsarr_out )
{
  if ( bitsarr_len < sub_sz )
    {
      cout << "Error! Input kmer too short." << bitsarr_len << " " << sub_sz << endl;
      return;
    }

  int arr_sz_in = bitsarr_len / 32 + 1;
  int rem = bitsarr_len % 32;

  if ( rem == 0 )
    {
      arr_sz_in--;
    }

  int arr_sz_out = sub_sz / 32 + 1;

  if ( sub_sz % 32 == 0 )
    {
      arr_sz_out--;
    }

  uint64_t temp_arr[10];
  memset ( temp_arr, 0, sizeof ( temp_arr ) );
  memset ( bitsarr_out, 0, sizeof ( uint64_t ) *arr_sz_out );
  int rem2 = ( 32 - rem + begin_pos ) % 32;
  int block_beg = ( 32 - rem + begin_pos ) / 32;

  if ( rem == 0 )
    {
      block_beg--;
    }

  int rem3 = ( 32 - rem + begin_pos + sub_sz ) % 32;
  int block_end = ( 32 - rem + begin_pos + sub_sz ) / 32;

  if ( rem3 != 0 )
    {
      rem3 = 32 - rem3;
    }
  else
    {
      block_end--;
    }

  if ( rem == 0 )
    {
      block_end--;
    }

  int orig_sz = ( block_end - block_beg + 1 );
  memcpy ( temp_arr, &bitsarr_in[block_beg], orig_sz * sizeof ( uint64_t ) );
  L_shift_NC ( temp_arr, rem2 * 2, orig_sz );
  R_shift_NC ( temp_arr, ( rem2 + rem3 ) % 32 * 2, arr_sz_out );
  memcpy ( bitsarr_out, temp_arr, arr_sz_out * sizeof ( uint64_t ) );
}

inline void L_shift_NC ( uint64_t *bitsarr, int shift_sz, int arr_sz )
{
  uint64_t temp_arr[100];

  for ( int i = 0; i < arr_sz; ++i )
    {
      temp_arr[i] = 0;
    }

  int jmp = shift_sz / 64;
  int offset = shift_sz % 64;

  for ( int i = 0; i < arr_sz; ++i )
    {
      if ( i + jmp + 1 < arr_sz )
        {
          uint64_t tt = 0;

          if ( offset == 0 )
            {
              tt = 0;
            }
          else
            {
              tt = ( bitsarr[i + jmp + 1] >> ( 64 - offset ) );
            }

          temp_arr[i] = ( ( bitsarr[i + jmp] << offset ) | tt );
        }

      if ( i + jmp + 1 == arr_sz )
        {
          temp_arr[i] = bitsarr[i + jmp] << offset;
        }

      if ( i + jmp + 1 > arr_sz )
        {
          temp_arr[i] = 0;
        }
    }

  for ( int i = 0; i < arr_sz; ++i )
    {
      bitsarr[i] = temp_arr[i];
    }
}



inline void R_shift_NC ( uint64_t *bitsarr, int shift_sz, int arr_sz )
{
  uint64_t temp_arr[100];

  for ( int i = 0; i < arr_sz; ++i )
    {
      temp_arr[i] = 0;
    }

  int jmp = shift_sz / 64;
  int offset = shift_sz % 64;

  if ( offset == 0 ) //to fix the move 64bit bug
    {
      for ( int i = arr_sz - 1; i >= 0; --i )
        {
          if ( i - jmp > 0 )
            {
              temp_arr[i] = bitsarr[i - jmp];
            }

          if ( i - jmp == 0 )
            {
              temp_arr[i] = bitsarr[i - jmp];
            }

          if ( i - jmp < 0 )
            {
              temp_arr[i] = 0;
            }
        }
    }
  else
    {
      for ( int i = arr_sz - 1; i >= 0; --i )
        {
          if ( i - jmp > 0 )
            {
              temp_arr[i] = ( bitsarr[i - jmp] >> offset ) | ( bitsarr[i - jmp - 1] << ( 64 - offset ) );
            }

          if ( i - jmp == 0 )
            {
              temp_arr[i] = ( bitsarr[i - jmp] >> offset );
            }

          if ( i - jmp < 0 )
            {
              temp_arr[i] = 0;
            }
        }
    }

  for ( int i = 0; i < arr_sz; ++i )
    {
      bitsarr[i] = temp_arr[i];
    }
}


inline int uint64_t_cmp ( uint64_t *A, uint64_t *B, int Kmer_arr_sz )
{
  int flag = 0;

  for ( int jj = 0; jj < Kmer_arr_sz; ++jj )
    {
      if ( A[jj] > B[jj] )
        {
          flag = 1;
          break;
        }

      if ( A[jj] < B[jj] )
        {
          flag = -1;
          break;
        }

      if ( A[jj] == B[jj] )
        {
          continue;
        }
    }

  return flag;
}


//for 63mer
inline uint64_t *get_rev_comp_seq_arr ( uint64_t *seq_arr, int seq_size, int arr_sz )
{
  if ( seq_size < 32 && arr_sz == 2 )
    {
      seq_arr[1] = get_rev_comp_seq ( seq_arr[1], seq_size );

      if ( seq_arr[0] != 0 )
        {
          fprintf ( stderr, "ERROR: in get_rev_comp_seq_arr \n" );
          exit ( -1 );
        }

      return seq_arr;
    }

  int tot_bits = arr_sz * 64;

  for ( int i = 0; i < arr_sz; ++i )
    {
      seq_arr[i] = ~seq_arr[i];
      seq_arr[i] = ( ( seq_arr[i] & 0x3333333333333333 ) << 2 ) | ( ( seq_arr[i] & 0xCCCCCCCCCCCCCCCC ) >> 2 );
      seq_arr[i] = ( ( seq_arr[i] & 0x0F0F0F0F0F0F0F0F ) << 4 ) | ( ( seq_arr[i] & 0xF0F0F0F0F0F0F0F0 ) >> 4 );
      seq_arr[i] = ( ( seq_arr[i] & 0x00FF00FF00FF00FF ) << 8 ) | ( ( seq_arr[i] & 0xFF00FF00FF00FF00 ) >> 8 );
      seq_arr[i] = ( ( seq_arr[i] & 0x0000FFFF0000FFFF ) << 16 ) | ( ( seq_arr[i] & 0xFFFF0000FFFF0000 ) >> 16 );
      seq_arr[i] = ( ( seq_arr[i] & 0x00000000FFFFFFFF ) << 32 ) | ( ( seq_arr[i] & 0xFFFFFFFF00000000 ) >> 32 );
    }

  int j = 0, k = arr_sz - 1;

  for ( ; j < k; ++j, --k )
    {
      uint64_t temp;
      temp = seq_arr[j];
      seq_arr[j] = seq_arr[k];
      seq_arr[k] = temp;
    }

  R_shift_NC ( seq_arr, tot_bits - ( seq_size * 2 ), arr_sz );
  return seq_arr;
}


inline uint64_t get_rev_comp_seq ( uint64_t seq, int seq_size )
{
  seq = ~seq;
  seq = ( ( seq & 0x3333333333333333 ) << 2 ) | ( ( seq & 0xCCCCCCCCCCCCCCCC ) >> 2 );
  seq = ( ( seq & 0x0F0F0F0F0F0F0F0F ) << 4 ) | ( ( seq & 0xF0F0F0F0F0F0F0F0 ) >> 4 );
  seq = ( ( seq & 0x00FF00FF00FF00FF ) << 8 ) | ( ( seq & 0xFF00FF00FF00FF00 ) >> 8 );
  seq = ( ( seq & 0x0000FFFF0000FFFF ) << 16 ) | ( ( seq & 0xFFFF0000FFFF0000 ) >> 16 );
  seq = ( ( seq & 0x00000000FFFFFFFF ) << 32 ) | ( ( seq & 0xFFFFFFFF00000000 ) >> 32 );
  return seq >> ( 64 - ( seq_size * 2 ) );
}


//for 64bit platform
inline uint64_t MurmurHash64A ( const void *key, int len, unsigned int seed )
{
  const uint64_t m = 0xc6a4a7935bd1e995;
  const int r = 47;
  uint64_t h = seed ^ ( len * m );
  const uint64_t *data = ( const uint64_t * ) key;
  const uint64_t *end = data + ( len / 8 );

  while ( data != end )
    {
      uint64_t k = *data++;
      k *= m;
      k ^= k >> r;
      k *= m;
      h ^= k;
      h *= m;
    }

  const unsigned char *data2 = ( const unsigned char * ) data;

  switch ( len & 7 )
    {
    case 7:
      h ^= uint64_t ( data2[6] ) << 48;

    case 6:
      h ^= uint64_t ( data2[5] ) << 40;

    case 5:
      h ^= uint64_t ( data2[4] ) << 32;

    case 4:
      h ^= uint64_t ( data2[3] ) << 24;

    case 3:
      h ^= uint64_t ( data2[2] ) << 16;

    case 2:
      h ^= uint64_t ( data2[1] ) << 8;

    case 1:
      h ^= uint64_t ( data2[0] );
      h *= m;
    };

  h ^= h >> r;

  h *= m;

  h ^= h >> r;

  return h;
}


//for 32bit platform
inline uint64_t MurmurHash64B ( const void *key, int len, unsigned int seed )
{
  const unsigned int m = 0x5bd1e995;
  const int r = 24;
  unsigned int h1 = seed ^ len;
  unsigned int h2 = 0;
  const unsigned int *data = ( const unsigned int * ) key;

  while ( len >= 8 )
    {
      unsigned int k1 = *data++;
      k1 *= m;
      k1 ^= k1 >> r;
      k1 *= m;
      h1 *= m;
      h1 ^= k1;
      len -= 4;
      unsigned int k2 = *data++;
      k2 *= m;
      k2 ^= k2 >> r;
      k2 *= m;
      h2 *= m;
      h2 ^= k2;
      len -= 4;
    }

  if ( len >= 4 )
    {
      unsigned int k1 = *data++;
      k1 *= m;
      k1 ^= k1 >> r;
      k1 *= m;
      h1 *= m;
      h1 ^= k1;
      len -= 4;
    }

  switch ( len )
    {
    case 3:
      h2 ^= ( ( unsigned char * ) data ) [2] << 16;

    case 2:
      h2 ^= ( ( unsigned char * ) data ) [1] << 8;

    case 1:
      h2 ^= ( ( unsigned char * ) data ) [0];
      h2 *= m;
    };

  h1 ^= h2 >> 18;

  h1 *= m;

  h2 ^= h1 >> 22;

  h2 *= m;

  h1 ^= h2 >> 17;

  h1 *= m;

  h2 ^= h1 >> 19;

  h2 *= m;

  uint64_t h = h1;

  h = ( h << 32 ) | h2;

  return h;
}

#endif


/*
 * convert_soapdenovo.cpp
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
#include "global.h"
#include "stdinc.h"
#include "core.h"


#include "convert_soapdenovo.h"
#include "sparse_kmer.h"
#include "zlib.h"

// can be optimized ...
/*************************************************
Function:
    convert_kmer
Description:
    convert sparse kmer to soapdenovo format ..
    A 0           A 0
    C 1  -->       C 1
    T 3             T 2
    G 2        G 3
Input:
    1. sparse_kmer:     sparse-kmer
    2. K_size:      kmer size
Output:
    1. sparse_kmer      soapdenovo format kmer
Return:
    None.
*************************************************/
void convert_kmer ( kmer_t2 * sparse_kmer, int K_size )
{
	uint64_t tmp, tmp_res;
	int i, j, index, arr_sz, base;
	index = K_size / 32;
	arr_sz = sizeof ( kmer_t2 ) / sizeof ( uint64_t );

	for ( i = 0; i <= index; i++ )
	{
		tmp = ( sparse_kmer->kmer ) [arr_sz - 1 - i];
		tmp_res = 0;

		for ( j = 0; j < 32; j++ )
		{
			base = tmp & 3;

			switch ( base )
			{
				case 0 :
					break;
				case 1 :
					tmp_res |= ( 1LLU << 2 * j );
					break;
				case 2:
					tmp_res |= ( 3LLU << 2 * j );
					break;
				case 3:
					tmp_res |= ( 2LLU << 2 * j );
					break;
			}

			tmp = tmp >> 2;
		}

		( sparse_kmer->kmer ) [arr_sz - 1 - i] = tmp_res;
	}

	/*
	uint64_t high=0,low=0;
	int chr;
	if(K_size>=33){
	    for(int i=0;i<K_size-32;++i){
	        chr = sparse_kmer[0] & 3;
	        switch(chr){
	            case 0:
	                break;
	            case 1:
	                high|=(1LLU << 2 * i);
	                break;
	            case 2:
	                high|=(3LLU << 2 * i);
	                break;
	            case 3:
	                high|=(2LLU << 2 * i);
	                break;
	        }

	        sparse_kmer[0] = sparse_kmer[0]>>2;
	    }

	    for(int i=0;i<32;++i){
	        chr = sparse_kmer[1] & 3;
	        switch(chr){
	            case 0:
	                break;
	            case 1:
	                low|=(1LLU << 2 * i);
	                break;
	            case 2:
	                low|=(3LLU << 2 * i);
	                break;
	            case 3:
	                low|=(2LLU << 2 * i);
	                break;
	        }

	        sparse_kmer[1] = sparse_kmer[1]>>2;
	    }
	}else{
	    for(int i=0;i<K_size;++i){
	        chr = sparse_kmer[1] & 3;
	        switch(chr){
	            case 0:
	                break;
	            case 1:
	                low|=(1LLU << 2 * i);
	                break;
	            case 2:
	                low|=(3LLU << 2 * i);
	                break;
	            case 3:
	                low|=(2LLU << 2 * i);
	                break;
	        }

	        sparse_kmer[1] = sparse_kmer[1]>>2;
	    }
	}
	sparse_kmer[0] = high;
	sparse_kmer[1] = low;
	*/
}


/*************************************************
Function:
     fastReverseComp
Description:
    fastReverseComp soapdenovo format kmer
Input:
    1. kmer2:       soapdenovo fomat kmer
    2. seq_size:        kmer size
Output:
    1. kmer2:       reversed compliment kmer
Return:
    None.
*************************************************/
static  void fastReverseComp ( kmer_t2 * kmer2, int seq_size )
{
	int arr_sz;
	uint64_t * seq_arr;
	arr_sz = sizeof ( kmer_t2 ) / sizeof ( uint64_t ); //= 2 or 4
	seq_arr = kmer2->kmer;
	int tot_bits = arr_sz * 64;

	for ( int i = 0; i < arr_sz; ++i )
	{
		seq_arr[i] ^= 0xAAAAAAAAAAAAAAAALLU;
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
}

struct classcomp
{
	bool operator() ( const kmer_t2 & t1, const kmer_t2 & t2 ) const
	{
		int Kmer_arr_sz = sizeof ( kmer_t2 ) / sizeof ( uint64_t );

		for ( int jj = 0; jj < Kmer_arr_sz; ++jj )
		{
			if ( ( t1.kmer ) [jj] < ( t2.kmer ) [jj] )
			{
				return 1;
			}
			else if ( ( t1.kmer ) [jj] > ( t2.kmer ) [jj] )
			{
				return 0;
			}

			continue;
		}

		return 0;
		/* old
		if((t1.kmer)[0] < (t2.kmer)[0]){
		    return 1;
		}else if((t1.kmer)[0] == (t2.kmer)[0]){
		    return (t1.kmer)[1] < (t2.kmer)[1];
		}else{
		    return 0;
		}*/
	}
};



/*************************************************
Function:
    convert
Description:
    converts *.sparse.edge to *.edge.gz *.vertex
Input:
    1. sparse_edge_file:        sparse edge file
    2. K_size:      kmer size
    3. output_prefix        output prefix
Output:
    None.
Return:
    None.
*************************************************/
void convert ( char * sparse_edge_file, int K_size, char * output_prefix )
{
	if ( run_mode != 0 )
	{
		char temp[256];
		sprintf ( temp, "%s.preGraphBasic", output_prefix );
		FILE * fp = fopen ( temp, "r" );
		char line[1024];
		fgets ( line, 1024, fp );
		fgets ( line, 1024, fp );
		fgets ( line, 1024, fp );
		fclose ( fp );
		sscanf ( line, "%s %d %s %d", temp, &max_rd_len, temp, &min_rd_len );
	}

	FILE * fin, *fout2, *fout3;
	fin = fopen ( sparse_edge_file, "r" );
	gzFile fout;
	char temp[256];
	sprintf ( temp, "%s.edge.gz", output_prefix );
	fout = gzopen ( temp, "w" );
	//fout= fopen(temp, "w");//edge
	//write as gzip file
	sprintf ( temp, "%s.vertex", output_prefix );
	fout2 = fopen ( temp, "w" );
	sprintf ( temp, "%s.preGraphBasic", output_prefix );
	fout3 = fopen ( temp, "w" );

	if ( !fin || !fout || !fout2 || !fout3 )
	{
		fprintf ( stderr,  "can't open file %s\n", sparse_edge_file );
		exit ( 1 );
	}

	//cout << "right 0"<<endl;
	kmer_t2 from_kmer, to_kmer;
	size_t line_len, edge_len_left;
	int edge_len;
	int cvg;
	int bal_ed;//»ØÎÄÎª0
	char str[32];
	const int BUFF_LEN = 1024;
	char line[BUFF_LEN];
	int start = 0;
	int cutoff = 100;
	map<kmer_t2, int, classcomp> vertex_nodes;
	size_t edge_counter = 0, vertex_counter = 0;
	int j = 0;

	//cout << "right 1"<<endl;

	while ( fgets ( line, BUFF_LEN, fin ) != NULL )
	{
		//cout << "right 2"<<endl;
		if ( line[0] == '>' ) //get one edge length, from vertex, to vertex,cvg,bal
		{
			edge_counter++;
#ifdef _63MER_
			sscanf ( line + 7, "%d,%llx %llx,%llx %llx,cvg %d,%d", &edge_len,
			         & ( from_kmer.kmer ) [0], & ( from_kmer.kmer ) [1], & ( to_kmer.kmer ) [0], & ( to_kmer.kmer ) [1], &cvg, &bal_ed ); // from_kmer to_kmer is of no use here
#endif
#ifdef _127MER_
			sscanf ( line + 7, "%d,%llx %llx %llx %llx,%llx %llx %llx %llx,cvg %d,%d",
			         &edge_len, & ( from_kmer.kmer ) [0], & ( from_kmer.kmer ) [1], & ( from_kmer.kmer ) [2], & ( from_kmer.kmer ) [3],
			         & ( to_kmer.kmer ) [0], & ( to_kmer.kmer ) [1], & ( to_kmer.kmer ) [2], & ( to_kmer.kmer ) [3], &cvg, &bal_ed ); // from_kmer to_kmer is of no use here
#endif

			if ( edge_len == 1 )
			{
				cvg = 0;
			}
			else
			{
				cvg *= 10;
			}

			convert_kmer ( &from_kmer, K_size );
			convert_kmer ( &to_kmer, K_size );
#ifdef _63MER_
			gzprintf ( fout, ">length %d,%llx %llx,%llx %llx,cvg %d,%d\n", edge_len,
			           ( from_kmer.kmer ) [0], ( from_kmer.kmer ) [1], ( to_kmer.kmer ) [0], ( to_kmer.kmer ) [1], cvg, bal_ed );
#endif
#ifdef _127MER_
			gzprintf ( fout, ">length %d,%llx %llx %llx %llx,%llx %llx %llx %llx,cvg %d,%d\n", edge_len,
			           ( from_kmer.kmer ) [0], ( from_kmer.kmer ) [1], ( from_kmer.kmer ) [2], ( from_kmer.kmer ) [3],
			           ( to_kmer.kmer ) [0], ( to_kmer.kmer ) [1], ( to_kmer.kmer ) [2], ( to_kmer.kmer ) [3], cvg, bal_ed );
#endif

			if ( bal_ed ) { edge_counter++; }

			kmer_t2 f_kmer = from_kmer;
			fastReverseComp ( &f_kmer, K_size );

			if ( kmerCompare ( &f_kmer, &from_kmer ) < 0 )
			{
				from_kmer = f_kmer;
			}

			vertex_nodes[from_kmer]++;
			f_kmer = to_kmer;
			fastReverseComp ( &f_kmer, K_size );

			if ( kmerCompare ( &f_kmer, &to_kmer ) < 0 )
			{
				to_kmer = f_kmer;
			}

			vertex_nodes[to_kmer]++;
			start = 1;
			j = 0;
		}
		else
		{
			//print the sequence
			if ( start == 1 )
			{
				//skip the first kmer
				int len = strlen ( line );

				if ( line[len - 1] == '\n' )
				{
					line[len - 1] == '\0';
					len --;
				}

				for ( int i = K_size; i < len; i++ )
				{
					j++;
					gzprintf ( fout, "%c", line[i] );

					if ( j % 100 == 0 )
					{
						gzprintf ( fout, "\n" );
					}
				}

				edge_len -= ( len - K_size );

				if ( edge_len == 0 && j % 100 != 0 )
				{
					gzprintf ( fout, "\n" );
				}

				start = 2;
			}
			else  //start = 2
			{
				if ( line[0] == '\n' ) { continue; }

				int len = strlen ( line );

				if ( line[len - 1] == '\n' )
				{
					line[len - 1] == '\0';
					len --;
				}

				for ( int i = 0; i < len; i++ )
				{
					j++;
					gzprintf ( fout, "%c", line[i] );

					if ( j % 100 == 0 )
					{
						gzprintf ( fout, "\n" );
					}
				}

				edge_len -= len;

				if ( edge_len == 0 && j % 100 != 0 )
				{
					gzprintf ( fout, "\n" );
				}
			}
		}
	}

	//fprintf(stderr,"size of map: %llu\n",vertex_nodes.size());
	map<kmer_t2, int>::iterator it;

	for ( it = vertex_nodes.begin(); it != vertex_nodes.end(); ++it )
	{
		vertex_counter++;
#ifdef _63MER_
		fprintf ( fout2, "%llx %llx ", ( ( *it ).first.kmer ) [0], ( ( *it ).first.kmer ) [1] );
#endif
#ifdef _127MER_
		fprintf ( fout2, "%llx %llx %llx %llx ", ( ( *it ).first.kmer ) [0], ( ( *it ).first.kmer ) [1],
		          ( ( *it ).first.kmer ) [2], ( ( *it ).first.kmer ) [3] );
#endif

		if ( vertex_counter % 8 == 0 ) { fprintf ( fout2, "\n" ); }
	}

	fprintf ( fout3, "VERTEX %lu K %d\n", vertex_counter, K_size );
	fprintf ( fout3, "EDGEs %lu\n", edge_counter );
	fprintf ( stderr, "%llu edges and %llu vertexes constructed.\n", edge_counter, vertex_counter );
	fprintf ( fout3, "MaxReadLen %d MinReadLen %d MaxNameLen 256\n", max_rd_len, min_rd_len );
	fclose ( fin );
	gzclose ( fout );
	fclose ( fout2 );
	fclose ( fout3 );
}



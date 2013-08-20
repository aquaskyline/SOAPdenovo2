/*
 * readseq1by1.c
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

#include <math.h>
#include "bam.h"
#include "faidx.h"
#include "knetfile.h"
#include "sam_view.h"
#include "xcurses.h"
#include "zlib.h"
#include "bgzf.h"
#include "glf.h"
#include "kstring.h"
#include "razf.h"
#include "sam_header.h"
#include "zconf.h"

static char src_rc_seq[1024];
static int state = -3;
static int readstate = 0;
boolean change = 1;
//static boolean multi=0;
//atic long curr_pos;

void readseq1by1 ( char * src_seq, char * src_name, int * len_seq, FILE * fp, long long num_seq )
{
	int i, k, n, strL;
	char c;
	int lineLen = gLineLen;
	char tmpStr[lineLen];
	char * str;

	if ( gStr == NULL )
	{
		str = tmpStr;
	}
	else
	{
		str = gStr;
		lineLen = maxReadLen + 1;
	}

	n = 0;
	k = num_seq;

	while ( fgets ( str, lineLen, fp ) )
	{
		if ( str[0] == '#' )
		{
			continue;
		}

		if ( str[0] == '>' )
		{
			/*
			   if(k >= 0)  {  // if this isn't the first '>' in the file
			   *len_seq = n;
			   }
			 */
			/*
			if(multi)
			{
			curr_pos=ftell(fp);
			}
			*/
			*len_seq = n;
			n = 0;
			sscanf ( &str[1], "%s", src_name );
			return;
		}
		else
		{
			strL = strlen ( str );

			if ( strL + n > maxReadLen )
			{
				strL = maxReadLen - n;
			}

			for ( i = 0; i < strL; i++ )
			{
				if ( str[i] >= 'a' && str[i] <= 'z' )
				{
					c = base2int ( str[i] - 'a' + 'A' );
					src_seq[n++] = c;
				}
				else if ( str[i] >= 'A' && str[i] <= 'Z' )
				{
					c = base2int ( str[i] );
					src_seq[n++] = c;
					// after pre-process all the symbles would be a,g,c,t,n in lower or upper case.
				}
				else if ( str[i] == '.' )
				{
					c = base2int ( 'A' );
					src_seq[n++] = c;
				}   // after pre-process all the symbles would be a,g,c,t,n in lower or upper case.
			}

			//printf("%d: %d\n",k,n);
		}
	}

	if ( k >= 0 )
	{
		*len_seq = n;
		return;
	}

	*len_seq = 0;
}

void readseqInBuf ( char * src_seq, char * src_name, int * len_seq, char * buf, int * start, int offset )
{
	int i, n, strL, m, p;
	char c;
	int lineLen = gLineLen;
	char tmpStr[lineLen];
	char * str;

	if ( gStr == NULL )
	{
		str = tmpStr;
	}
	else
	{
		str = gStr;
		lineLen = maxReadLen + 1;
	}

	n = 0;

	for ( m = *start; m < offset; m++ )
	{
		if ( buf[m] == '>' )
		{
			p = m;
		}
		else if ( buf[m] == '\n' && buf[p] == '>' )
		{
			memcpy ( src_name, &buf[p + 1], m - p - 1 ); //get name
			p = m;
		}
		else if ( buf[m] == '\n' && buf[p] == '\n' )
		{
			memcpy ( str, &buf[p + 1], m - p - 1 ); //get seq
			//p = m;
			str[m - p - 1] = '\0';
			*start = m + 1;
			strL = strlen ( str );

			if ( strL + n > maxReadLen )
				{ strL = maxReadLen - n; }

			for ( i = 0; i < strL; i++ )
			{
				if ( str[i] >= 'a' && str[i] <= 'z' )
				{
					c = base2int ( str[i] - 'a' + 'A' );
					src_seq[n++] = c;
				}
				else if ( str[i] >= 'A' && str[i] <= 'Z' )
				{
					c = base2int ( str[i] );
					src_seq[n++] = c;
					// after pre-process all the symbles would be a,g,c,t,n in lower or upper case.
				}
				else if ( str[i] == '.' )
				{
					c = base2int ( 'A' );
					src_seq[n++] = c;
				}   // after pre-process all the symbles would be a,g,c,t,n in lower or upper case.
			}

			break;
			//printf("%d: %d\n",k,n);
		}
	}

	*len_seq = n;
	return;
}

/*************************************************
Function:
    readseqpar
Description:
    Parses a contig file  to get the max contig length, the min contig length, the max name length.
Input:
    1. fp:      the contig file
Output:
    1. max_leg:         the max contig length
    2. min_leg:         the min contig length
    3. max_name_leg:        the max name length
Return:
    The contig number.
*************************************************/
long long readseqpar ( int * max_leg, int * min_leg, int * max_name_leg, FILE * fp )
{
	int l, n;
	long long k;
	char str[5000], src_name[5000];
	n = 0;
	k = -1;

	while ( fgets ( str, 4950, fp ) )
	{
		if ( str[0] == '>' )
		{
			if ( k >= 0 )
			{
				if ( n > *max_leg )
				{
					*max_leg = n;
				}

				if ( n < *min_leg )
				{
					*min_leg = n;
				}
			}

			n = 0;
			k++;
			sscanf ( &str[1], "%s", src_name );

			if ( ( l = strlen ( src_name ) ) > *max_name_leg )
			{
				*max_name_leg = l;
			}
		}
		else
		{
			n += strlen ( str ) - 1;
		}
	}

	if ( n > *max_leg )
	{
		*max_leg = n;
	}

	if ( n < *min_leg )
	{
		*min_leg = n;
	}

	k++;
	return ( k );
}

void readseqfq ( char * src_seq, char * src_name, int * len_seq, char * buf, int * start, int offset )
{
	int i, n, strL, m, p = 0;
	char c;
	int lineLen = gLineLen;
	char tmpStr[lineLen];
	char * str;

	if ( gStr == NULL )
	{
		str = tmpStr;
	}
	else
	{
		str = gStr;
		lineLen = maxReadLen + 1;
	}

	n = 0;

	for ( m = *start; m < offset; m++ )
	{
		if ( buf[m] == '@' )
			{ p = m; }
		else if ( buf[m] == '\n' && buf[p] == '@' )
		{
			memcpy ( src_name, &buf[p + 1], m - p - 1 );
			p = m;
		}
		else if ( buf[m] == '\n' && buf[p] == '\n' && m > p )
		{
			//p = m;
			memcpy ( str, &buf[p + 1], m - p - 1 );
			str[m - p - 1] = '\0';
			strL = strlen ( str );

			if ( strL + n > maxReadLen )
				{ strL = maxReadLen - n; }

			for ( i = 0; i < strL; i++ )
			{
				if ( str[i] >= 'a' && str[i] <= 'z' )
				{
					c = base2int ( str[i] - 'a' + 'A' );
					src_seq[n++] = c;
				}
				else if ( str[i] >= 'A' && str[i] <= 'Z' )
				{
					c = base2int ( str[i] );
					src_seq[n++] = c;
					// after pre-process all the symbles would be a,g,c,t,n in lower or upper case.
				}
				else if ( str[i] == '.' )
				{
					c = base2int ( 'A' );
					src_seq[n++] = c;
				}
			}

			for ( m++; buf[m] != '\n'; m++ ) { ; }

			*start = ( m + 2 + strlen ( str ) );
			break;
		}

		/*  else if (buf[m] == '\n' && buf[p] == '+' && buf[m - 1] != '+')
		    {
		        *start = m + 1;
		        break;
		    }*/
	}

	*len_seq = n;
	return;
}

void read1seqfq ( char * src_seq, char * src_name, int * len_seq, FILE * fp )
{
	int i, n, strL;
	char c;
	int lineLen = gLineLen;
	char tmpStr[lineLen];
	char * str;

	if ( gStr == NULL )
	{
		str = tmpStr;
	}
	else
	{
		str = gStr;
		lineLen = maxReadLen + 1;
	}

	boolean flag = 0;

	while ( fgets ( str, lineLen, fp ) )
	{
		if ( str[0] == '@' )
		{
			flag = 1;
			sscanf ( &str[1], "%s", src_name );
			break;
		}
	}

	if ( !flag )    //last time reading fq file get this
	{
		*len_seq = 0;
		return;
	}

	/*
	if(multi)
	{
	    curr_pos=ftell(fp);
	}
	*/
	n = 0;

	while ( fgets ( str, lineLen, fp ) )
	{
		if ( str[0] == '+' )
		{
			fgets ( str, lineLen, fp ); // pass quality value line
			*len_seq = n;
			return;
		}
		else
		{
			strL = strlen ( str );

			if ( strL + n > maxReadLen )
			{
				strL = maxReadLen - n;
			}

			for ( i = 0; i < strL; i++ )
			{
				if ( str[i] >= 'a' && str[i] <= 'z' )
				{
					c = base2int ( str[i] - 'a' + 'A' );
					src_seq[n++] = c;
				}
				else if ( str[i] >= 'A' && str[i] <= 'Z' )
				{
					c = base2int ( str[i] );
					src_seq[n++] = c;
					// after pre-process all the symbles would be a,g,c,t,n in lower or upper case.
				}
				else if ( str[i] == '.' )
				{
					c = base2int ( 'A' );
					src_seq[n++] = c;
				}   // after pre-process all the symbles would be a,g,c,t,n in lower or upper case.
			}
		}
	}

	*len_seq = n;
	return;
}

void read1seqbam ( char * src_seq, char * src_name, int * len_seq, samfile_t * in, int * type, int asm_flag ) //read one sequence from bam file
{
	bam1_t * b = bam_init1 ();
	char c;
	char * line1 = NULL;
	int n = 0;
	int len;
	int i, j;
	char * seq1;
	unsigned int flag1 = 0;
	*type = 0;
	readstate = 0;
	boolean isGood = true;

	//  if ((readstate = samread (in, b)) >= 0)
	while ( ( readstate = samread ( in, b ) ) >= 0 )
	{
		if ( !__g_skip_aln ( in->header, b ) )
		{
			line1 = bam_format1_core ( in->header, b, in->type >> 2 & 3 );
		}

		seq1 = strtok ( line1, "\t" );

		for ( i = 0; i < 10; i++ )
		{
			if ( i == 0 )
			{
				sscanf ( seq1, "%s", src_name );
			}
			else if ( i == 1 )
			{
				flag1 = atoi ( seq1 );

				if ( flag1 & 0x0200 ) //whether it's good or not
				{
					if ( asm_flag == 1 )
					{
						isGood = false;
						break;
					}

					switch ( state )
					{
						case -3:
							state = -2;
							break;
						case -2:
							state = 0;
							break;
						case -1:
							state = 2;
							break;
						default:
							state = -3;
					}
				}
				else
				{
					isGood = true;

					switch ( state )
					{
						case -3:
							state = -1;
							break;
						case -2:
							state = 1;
							break;
						case -1:
							state = 3;
							break;
						default:
							state = -3;
					}
				}
			}
			else if ( i == 9 )  //the sequence
			{
				len = strlen ( seq1 );

				if ( len + n > maxReadLen )
					{ len = maxReadLen - n; }

				for ( j = 0; j < len; j++ )
				{
					if ( seq1[j] >= 'a' && seq1[j] <= 'z' )
					{
						c = base2int ( seq1[j] - 'a' + 'A' );
						src_seq[n++] = c;
					}
					else if ( seq1[j] >= 'A' && seq1[j] <= 'Z' )
					{
						c = base2int ( seq1[j] );
						src_seq[n++] = c;
						// after pre-process all the symbles would be a,g,c,t,n in lower or upper case.
					}
					else if ( seq1[j] == '.' )
					{
						c = base2int ( 'A' );
						src_seq[n++] = c;
					}   // after pre-process all the symbles would be a,g,c,t,n in lower or upper case.
				}

				if ( 3 == state )
				{
					state = -3;
				}
				else
				{
					if ( 0 == state || 1 == state || 2 == state )
					{
						state = -3;
						*type = -1;
					}
				}
			}

			seq1 = strtok ( NULL, "\t" );
		}

		if ( isGood )
			{ break; }
	}

	if ( readstate < 0 )
		{ state = -3; }

	free ( line1 );
	bam_destroy1 ( b );
	*len_seq = n;
}

// find the next file to open in libs
int nextValidIndex ( int libNo, boolean pair, unsigned char asm_ctg )
{
	int i = libNo;

	while ( i < num_libs )
	{
		if ( asm_ctg == 1 && ( lib_array[i].asm_flag != 1 && lib_array[i].asm_flag != 3 ) )
		{
			i++;
			continue;
		}
		else if ( asm_ctg == 0 && ( lib_array[i].asm_flag != 2 && lib_array[i].asm_flag != 3 ) )
		{
			i++;
			continue;
		}
		else if ( asm_ctg > 1 && lib_array[i].asm_flag != asm_ctg ) // reads for other purpose
		{
			i++;
			continue;
		}

		if ( lib_array[i].curr_type == 1 && lib_array[i].curr_index < lib_array[i].num_a1_file )
		{
			return i;
		}

		if ( lib_array[i].curr_type == 2 && lib_array[i].curr_index < lib_array[i].num_q1_file )
		{
			return i;
		}

		if ( lib_array[i].curr_type == 3 && lib_array[i].curr_index < lib_array[i].num_p_file )
		{
			return i;
		}

		if ( lib_array[i].curr_type == 4 && lib_array[i].curr_index < lib_array[i].num_b_file )
		{
			return i;
		}

		if ( pair )
		{
			if ( lib_array[i].curr_type < 4 )
			{
				lib_array[i].curr_type++;
				lib_array[i].curr_index = 0;
			}
			else
			{
				i++;
			}

			continue;
		}

		if ( lib_array[i].curr_type == 5 && lib_array[i].curr_index < lib_array[i].num_s_a_file )
		{
			return i;
		}

		if ( lib_array[i].curr_type == 6 && lib_array[i].curr_index < lib_array[i].num_s_q_file )
		{
			return i;
		}

		if ( lib_array[i].curr_type < 6 )
		{
			lib_array[i].curr_type++;
			lib_array[i].curr_index = 0;
		}
		else
		{
			i++;
		}
	}           //for each lib

	return i;
}

static FILE * openFile4read ( char * fname )
{
	FILE * fp;
	int i, j = 0;
	char str_sub[4];

	for ( i = strlen ( fname ) - 1; i >= 0; i-- )
	{
		if ( fname[i] == ' ' )  { j++; }
		else { break; }
	}

	for ( i = 0; i < 3; i++ ) { str_sub[i] = fname[strlen ( fname ) - 3 - j + i]; }

	str_sub[3] = '\0';

	if ( strlen ( fname ) > 3 && strcmp ( str_sub, ".gz" ) == 0 )
	{
		char * cmd = ( char * ) ckalloc ( ( strlen ( fname ) + 20 ) * sizeof ( char ) );
		sprintf ( cmd, "gzip -dc %s", fname );
		fp = popen ( cmd, "r" );
		free ( cmd );
		return fp;
	}
	else
	{
		fname[strlen ( fname ) - j] = '\0';
		return ckopen ( fname, "r" );
	}
}

static samfile_t * openFile4readb ( char * fname ) //open file to read bam file
{
	samfile_t * in;
	char * fn_list = 0;

	if ( ( in = ( samfile_t * ) samopen ( fname, "rb", fn_list ) ) == 0 )
	{
		fprintf ( stderr, "Cannot open %s. Now exit to system...\n", fname );
		exit ( -1 );
	}

	if ( in->header == 0 )
	{
		fprintf ( stderr, "Cannot read the header.\n" );
		exit ( -1 );
	}

	return ( in );
}

void openFileInLib ( int libNo )
{
	int i = libNo;

	if ( lib_array[i].curr_type == 1 )
	{
		fprintf ( stderr, "Import reads from file:\n %s\n", lib_array[i].a1_fname[lib_array[i].curr_index] );
		fprintf ( stderr, "Import reads from file:\n %s\n", lib_array[i].a2_fname[lib_array[i].curr_index] );
		lib_array[i].fp1 = openFile4read ( lib_array[i].a1_fname[lib_array[i].curr_index] );
		lib_array[i].fp2 = openFile4read ( lib_array[i].a2_fname[lib_array[i].curr_index] );
		lib_array[i].curr_index++;
		lib_array[i].paired = 1;
	}
	else if ( lib_array[i].curr_type == 2 )
	{
		fprintf ( stderr, "Import reads from file:\n %s\n", lib_array[i].q1_fname[lib_array[i].curr_index] );
		fprintf ( stderr, "Import reads from file:\n %s\n", lib_array[i].q2_fname[lib_array[i].curr_index] );
		lib_array[i].fp1 = openFile4read ( lib_array[i].q1_fname[lib_array[i].curr_index] );
		lib_array[i].fp2 = openFile4read ( lib_array[i].q2_fname[lib_array[i].curr_index] );
		lib_array[i].curr_index++;
		lib_array[i].paired = 1;
	}
	else if ( lib_array[i].curr_type == 3 )
	{
		fprintf ( stderr, "Import reads from file:\n %s\n", lib_array[i].p_fname[lib_array[i].curr_index] );
		lib_array[i].fp1 = openFile4read ( lib_array[i].p_fname[lib_array[i].curr_index] );
		lib_array[i].curr_index++;
		lib_array[i].paired = 0;
	}
	else if ( lib_array[i].curr_type == 5 )
	{
		fprintf ( stderr, "Import reads from file:\n %s\n", lib_array[i].s_a_fname[lib_array[i].curr_index] );
		lib_array[i].fp1 = openFile4read ( lib_array[i].s_a_fname[lib_array[i].curr_index] );
		lib_array[i].curr_index++;
		lib_array[i].paired = 0;
	}
	else if ( lib_array[i].curr_type == 6 )
	{
		fprintf ( stderr, "Import reads from file:\n %s\n", lib_array[i].s_q_fname[lib_array[i].curr_index] );
		lib_array[i].fp1 = openFile4read ( lib_array[i].s_q_fname[lib_array[i].curr_index] );
		lib_array[i].curr_index++;
		lib_array[i].paired = 0;
	}
	else if ( lib_array[i].curr_type == 4 )
	{
		fprintf ( stderr, "Import reads from file:\n %s\n", lib_array[i].b_fname[lib_array[i].curr_index] );
		lib_array[i].fp3 = openFile4readb ( lib_array[i].b_fname[lib_array[i].curr_index] );
		lib_array[i].curr_index++;
		lib_array[i].paired = 0;
	}
}

static void reverse2k ( char * src_seq, int len_seq )
{
	if ( !len_seq )
	{
		return;
	}

	int i;
	reverseComplementSeq ( src_seq, len_seq, src_rc_seq );

	for ( i = 0; i < len_seq; i++ )
	{
		src_seq[i] = src_rc_seq[i];
	}
}

void closeFp1InLab ( int libNo )
{
	int ftype = lib_array[libNo].curr_type;
	int index = lib_array[libNo].curr_index - 1;
	char * fname;
	int i, j = 0;
	char str_sub[4];

	if ( ftype == 1 )
	{
		fname = lib_array[libNo].a1_fname[index];
	}
	else if ( ftype == 2 )
	{
		fname = lib_array[libNo].q1_fname[index];
	}
	else if ( ftype == 3 )
	{
		fname = lib_array[libNo].p_fname[index];
	}
	else if ( ftype == 5 )
	{
		fname = lib_array[libNo].s_a_fname[index];
	}
	else if ( ftype == 6 )
	{
		fname = lib_array[libNo].s_q_fname[index];
	}
	else if ( ftype == 4 )
	{
		fname = lib_array[libNo].b_fname[index];
	}
	else
	{
		return;
	}

	if ( ftype == 4 )
	{
		samclose ( lib_array[libNo].fp3 );  //close file3
	}
	else
	{
		for ( i = strlen ( fname ) - 1; i >= 0; i-- )
		{
			if ( fname[i] == ' ' )  { j++; }
			else { break; }
		}

		for ( i = 0; i < 3; i++ ) { str_sub[i] = fname[strlen ( fname ) - 3 - j + i]; }

		str_sub[3] = '\0';

		if ( strlen ( fname ) > 3 && strcmp ( str_sub, ".gz" ) == 0 )
		{
			pclose ( lib_array[libNo].fp1 );
		}
		else
		{
			fclose ( lib_array[libNo].fp1 );
		}
	}
}

void closeFp2InLab ( int libNo )
{
	int ftype = lib_array[libNo].curr_type;
	int index = lib_array[libNo].curr_index - 1;
	char * fname;
	int i, j = 0;
	char str_sub[4];

	if ( ftype == 1 )
	{
		fname = lib_array[libNo].a2_fname[index];
	}
	else if ( ftype == 2 )
	{
		fname = lib_array[libNo].q2_fname[index];
	}
	else
	{
		return;
	}

	for ( i = strlen ( fname ) - 1; i >= 0; i-- )
	{
		if ( fname[i] == ' ' )  { j++; }
		else { break; }
	}

	for ( i = 0; i < 3; i++ ) { str_sub[i] = fname[strlen ( fname ) - 3 - j + i]; }

	str_sub[3] = '\0';

	if ( strlen ( fname ) > 3 && strcmp ( str_sub, ".gz" ) == 0 )
	{
		pclose ( lib_array[libNo].fp2 );
	}
	else
	{
		fclose ( lib_array[libNo].fp2 );
	}
}

boolean readseqInLib ( char * src_seq, char * src_name, int * len_seq, char * buf, int * start, int offset, int i )
{
	if ( lib_array[i].curr_type == 1 )
	{
		if ( lib_array[i].paired == 1 )
		{
			readseqInBuf ( src_seq, src_name, len_seq, buf, start, offset );

			if ( lib_array[i].reverse )
				{ reverse2k ( src_seq, *len_seq ); }

			lib_array[i].paired = 2;
			n_solexa++;
			return 1;
		}
		else
		{
			readseqInBuf ( src_seq, src_name, len_seq, buf, start, offset );

			if ( lib_array[i].reverse )
				{ reverse2k ( src_seq, *len_seq ); }

			lib_array[i].paired = 1;
			n_solexa++;
			return 1;   //can't fail to read a read2
		}
	}

	if ( lib_array[i].curr_type == 2 )
	{
		if ( lib_array[i].paired == 1 )
		{
			readseqfq ( src_seq, src_name, len_seq, buf, start, offset );

			if ( lib_array[i].reverse )
				{ reverse2k ( src_seq, *len_seq ); }

			lib_array[i].paired = 2;
			n_solexa++;
			return 1;
		}
		else
		{
			readseqfq ( src_seq, src_name, len_seq, buf, start, offset );

			if ( lib_array[i].reverse )
				{ reverse2k ( src_seq, *len_seq ); }

			lib_array[i].paired = 1;
			n_solexa++;
			return 1;   //can't fail to read a read2
		}
	}

	if ( lib_array[i].curr_type == 6 )
		{ readseqfq ( src_seq, src_name, len_seq, buf, start, offset ); }
	else
	{
		readseqInBuf ( src_seq, src_name, len_seq, buf, start, offset );
	}

	/*
	   int t;
	   for(t=0;t<*len_seq;t++)
	   printf("%d",src_seq[t]);
	   printf("\n");
	 */
	if ( lib_array[i].reverse )
		{ reverse2k ( src_seq, *len_seq ); }

	n_solexa++;
	return 1;
}

boolean read1seqInLib ( char * src_seq, char * src_name, int * len_seq, int * libNo, boolean pair, unsigned char asm_ctg, int * type )
{
	int i = *libNo;
	int prevLib = i;

	if ( ( ( lib_array[i].curr_type != 4 ) && !lib_array[i].fp1 ) // file1 does not exist
	        || ( ( lib_array[i].curr_type == 4 ) && ( lib_array[i].fp3 == NULL ) ) //file3 does not exist
	        || ( ( lib_array[i].curr_type == 4 ) && readstate < 0 ) //file3 reaches end
	        || ( ( lib_array[i].curr_type == 1 ) && feof ( lib_array[i].fp1 ) && feof ( lib_array[i].fp2 ) ) //f1 && f2 reach end
		|| ( ( lib_array[i].curr_type == 2 ) && ( feof ( lib_array[i].fp1 ) || feof ( lib_array[i].fp2 ) )) //f1||f2 reaches end
	        || ( ( lib_array[i].curr_type != 1 && lib_array[i].curr_type != 2 ) && ( lib_array[i].curr_type != 4 ) && feof ( lib_array[i].fp1 ) ) ) // file1 reaches end and not type1 type2 and not type6
	{
		if ( lib_array[i].curr_type == 4 )
		{
			if ( lib_array[i].fp3 && readstate < 0 ) // file3 reaches end
				{ closeFp1InLab ( i ); }

			readstate = 0;
		}
		else
		{
			if ( lib_array[i].fp1 && feof ( lib_array[i].fp1 ) )
			{ 
				closeFp1InLab ( i );

				if ( lib_array[i].fp2 )
				{
					closeFp2InLab ( i );
				}
			}
			else if ( lib_array[i].fp2 && feof ( lib_array[i].fp2 ) )
			{
				closeFp2InLab ( i );

				if ( lib_array[i].fp1 )
				{
					closeFp1InLab ( i );
				}
			}
		}

		*libNo = nextValidIndex ( i, pair, asm_ctg );
		i = *libNo;

		if ( lib_array[i].rd_len_cutoff > 0 )
			{ maxReadLen = lib_array[i].rd_len_cutoff < maxReadLen4all ? lib_array[i].rd_len_cutoff : maxReadLen4all; }
		else
		{
			maxReadLen = maxReadLen4all;
		}

		//record insert size info
		if ( pair && i != prevLib )
		{
			if ( readNumBack < n_solexa )
			{
				pes[gradsCounter].PE_bound = n_solexa;
				pes[gradsCounter].rank = lib_array[prevLib].rank;
				pes[gradsCounter].pair_num_cut = lib_array[prevLib].pair_num_cut;
				pes[gradsCounter++].insertS = lib_array[prevLib].avg_ins;
				readNumBack = n_solexa;
			}
		}

		if ( i >= num_libs )
		{
			return 0;
		}

		openFileInLib ( i );

		if ( lib_array[i].curr_type == 1 )
		{
			readseq1by1 ( src_seq, src_name, len_seq, lib_array[i].fp1, -1 );
			readseq1by1 ( src_seq, src_name, len_seq, lib_array[i].fp2, -1 );
		}
		else if ( lib_array[i].curr_type == 3 || lib_array[i].curr_type == 5 )
		{
			readseq1by1 ( src_seq, src_name, len_seq, lib_array[i].fp1, -1 );
		}
	}

	if ( lib_array[i].curr_type == 1 )
	{
		if ( lib_array[i].paired == 1 )
		{
			readseq1by1 ( src_seq, src_name, len_seq, lib_array[i].fp1, 1 );

			if ( lib_array[i].reverse )
			{
				reverse2k ( src_seq, *len_seq );
			}

			lib_array[i].paired = 2;

			if ( *len_seq > 0 || !feof ( lib_array[i].fp1 ) )
			{
				n_solexa++;
				return 1;
			}
			else
			{
				return read1seqInLib ( src_seq, src_name, len_seq, libNo, pair, asm_ctg, type );
			}
		}
		else
		{
			readseq1by1 ( src_seq, src_name, len_seq, lib_array[i].fp2, 1 );

			if ( lib_array[i].reverse )
			{
				reverse2k ( src_seq, *len_seq );
			}

			lib_array[i].paired = 1;
			n_solexa++;
			return 1;   //can't fail to read a read2
		}
	}

	if ( lib_array[i].curr_type == 2 )
	{
		if ( lib_array[i].paired == 1 )
		{
			read1seqfq ( src_seq, src_name, len_seq, lib_array[i].fp1 );

			/*
			   if(*len_seq>0){
			   for(j=0;j<*len_seq;j++)
			   printf("%c",int2base(src_seq[j]));
			   printf("\n");
			   }
			 */
			if ( lib_array[i].reverse )
			{
				reverse2k ( src_seq, *len_seq );
			}

			lib_array[i].paired = 2;

			if ( *len_seq > 0 || !feof ( lib_array[i].fp1 ) )
			{
				n_solexa++;
				return 1;
			}
			else
			{
				return read1seqInLib ( src_seq, src_name, len_seq, libNo, pair, asm_ctg, type );
			}
		}
		else
		{
			read1seqfq ( src_seq, src_name, len_seq, lib_array[i].fp2 );

			if ( lib_array[i].reverse )
			{
				reverse2k ( src_seq, *len_seq );
			}

			lib_array[i].paired = 1;
			n_solexa++;
			return 1;   //can't fail to read a read2
		}
	}

	if ( lib_array[i].curr_type == 6 )
	{
		read1seqfq ( src_seq, src_name, len_seq, lib_array[i].fp1 );
	}
	else if ( lib_array[i].curr_type == 4 )
	{
		read1seqbam ( src_seq, src_name, len_seq, lib_array[i].fp3, type, lib_array[i].asm_flag );
	}
	else
	{
		readseq1by1 ( src_seq, src_name, len_seq, lib_array[i].fp1, 1 );
	}

	/*
	   int t;
	   for(t=0;t<*len_seq;t++)
	   printf("%d",src_seq[t]);
	   printf("\n");
	 */
	if ( lib_array[i].reverse )
	{
		reverse2k ( src_seq, *len_seq );
	}

	if ( lib_array[i].curr_type != 4 && ( *len_seq > 0 || !feof ( lib_array[i].fp1 ) ) )
	{
		n_solexa++;
		return 1;
	}
	else if ( lib_array[i].curr_type == 4 && ( *len_seq > 0 || readstate >= 0 ) )
	{
		n_solexa++;
		return 1;
	}
	else
	{
		return read1seqInLib ( src_seq, src_name, len_seq, libNo, pair, asm_ctg, type );
	}
}

boolean read1seqInLibBam ( char * src_seq, char * src_name, int * len_seq, int * libNo, boolean pair, unsigned char asm_ctg, int * type )
{
	int i = *libNo;

	if ( lib_array[i].fp3 == NULL )
	{
		fprintf ( stderr, "Empty file handle.\n" );
		return 0;
	}

	if ( lib_array[i].fp3 && readstate < 0 ) // file3 reaches end
	{
		closeFp1InLab ( i );
		readstate = 0;
		return 0;
	}

	read1seqbam ( src_seq, src_name, len_seq, lib_array[i].fp3, type, lib_array[i].asm_flag );

	if ( lib_array[i].reverse )
	{
		reverse2k ( src_seq, *len_seq );
	}

	if ( *len_seq > 0 || readstate >= 0 )
	{
		n_solexa++;
		return 1;
	}
	else
		{ return read1seqInLibBam ( src_seq, src_name, len_seq, libNo, pair, asm_ctg, type ); }
}

FILE * file = NULL;
boolean read1seqInLibpos ( char * src_seq, char * src_name, int * len_seq, // FILE *file,
                           int * file_No, int file_num, char ** fileName, int * fileType, int * maxLen, long * pos_seq )
{
	if ( *file_No < 0 ||  feof ( file ) )
	{
		( *file_No ) ++;

		if ( *file_No >= file_num )
			{ return 0; }

		maxReadLen = maxLen[*file_No];

		if ( file != NULL )
		{
			fclose ( file );
		}

		file = openFile4read ( fileName[*file_No] );

		if ( file != NULL )
			{ fprintf ( stderr, "Import reads from file:\n %s\n", fileName[*file_No] ); }

		if ( fileType[*file_No] == 1 || fileType[*file_No] == 3 )
		{
			readseq1by1 ( src_seq, src_name, len_seq, file, 1 );
			//          *pos_seq = curr_pos;
		}
	}

	//  multi=1;
	*len_seq = 0;

	if ( fileType[*file_No] == 1 )
	{
		//      *pos_seq = curr_pos;
		readseq1by1 ( src_seq, src_name, len_seq, file, 1 );
	}
	else if ( fileType[*file_No] == 2 )
	{
		read1seqfq ( src_seq, src_name, len_seq, file );
		//      *pos_seq = curr_pos;
	}
	else if ( fileType[*file_No] == 3 )
	{
		//      *pos_seq = curr_pos;
		readseq1by1 ( src_seq, src_name, len_seq, file, 1 );
	}
	else if ( fileType[*file_No] == 4 )
	{
		read1seqfq ( src_seq, src_name, len_seq, file );
		//      *pos_seq = curr_pos;
	}
	else
		{ return 0; }

	//  multi=0;
	if ( *len_seq > 0 || feof ( file ) )
	{
		n_solexa++;
		return 1;
	}
	else
		{ return read1seqInLibpos ( src_seq, src_name, len_seq, file_No, file_num, fileName, fileType, maxLen, pos_seq ); }
}


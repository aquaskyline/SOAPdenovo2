/*
 * check.c
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

#include <stdinc.h>
void * ckalloc ( unsigned long long amount );
FILE * ckopen ( char * name, char * mode );
FILE * ckopen ( char * name, char * mode )
{
	FILE * fp;

	if ( ( fp = fopen ( name, mode ) ) == NULL )
	{
		fprintf ( stderr, "Cannot open %s. Now exit to system...\n", name );
		exit ( -1 );
	}

	return ( fp );
}

static int GetFileSize ( FILE * fp )
{
	char c = fgetc ( fp );

	if ( c == EOF )
	{
		return 0;
	}

	return 1;
}

boolean check_file ( char * name )
{
	FILE * linkF;

	if ( strlen ( name ) > 3 && strcmp ( name + strlen ( name ) - 3, ".gz" ) == 0 )
	{
		char cmd[1000];
		sprintf ( cmd, "gzip -dc %s", name );
		linkF = popen ( cmd, "r" );

		if ( linkF )
		{
			if ( GetFileSize ( linkF ) != 0 )
			{
				pclose ( linkF );
				return 1;
			}

			pclose ( linkF );
		}

		return 0;
	}
	else
	{
		linkF = fopen ( name, "r" );

		if ( linkF )
		{
			if ( GetFileSize ( linkF ) != 0 )
			{
				fclose ( linkF );
				return 1;
			}

			fclose ( linkF );
		}

		return 0;
	}
}


/* ckalloc - allocate space; check for success */
void * ckalloc ( unsigned long long amount )
{
	void * p;

	if ( ( p = ( void * ) calloc ( 1, ( unsigned long long ) amount ) ) == NULL && amount != 0 )
	{
		fprintf ( stderr, "Ran out of memory while applying %lldbytes\n", amount );
		fprintf ( stderr, "There may be errors as follows:\n" );
		fprintf ( stderr, "1) Not enough memory.\n" );
		fprintf ( stderr, "2) The ARRAY may be overrode.\n" );
		fprintf ( stderr, "3) The wild pointers.\n" );
		exit ( -1 );
	}

	return ( p );
}


/* reallocate memory */
void * ckrealloc ( void * p, size_t new_size, size_t old_size )
{
	void * q;
	q = realloc ( ( void * ) p, new_size );

	if ( new_size == 0 || q != ( void * ) 0 )
	{
		return q;
	}

	/* manually reallocate space */
	q = ckalloc ( new_size );
	/* move old memory to new space */
	bcopy ( p, q, old_size );
	free ( p );
	return q;
}



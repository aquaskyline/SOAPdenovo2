/*
 * 31mer/output_pregraph.c
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
#include "newhash.h"
#include "kmerhash.h"
#include <extfunc.h>
#include <extvab.h>
static int outvCounter = 0;

//after this LINKFLAGFILTER in the Kmer is destroyed
static void output1vt ( kmer_t * node1, FILE * fp )
{
	if ( !node1 )
	{
		return;
	}

	if ( ! ( node1->linear ) && ! ( node1->deleted ) )
	{
		outvCounter++;
		print_kmer ( fp, node1->seq, ' ' );

		if ( outvCounter % 8 == 0 )
		{
			fprintf ( fp, "\n" );
		}
	}
}

void output_vertex ( char * outfile )
{
	char temp[256];
	FILE * fp;
	int i;
	kmer_t * node;
	KmerSet * set;
	sprintf ( temp, "%s.vertex", outfile );
	fp = ckopen ( temp, "w" );

	for ( i = 0; i < thrd_num; i++ )
	{
		set = KmerSets[i];
		set->iter_ptr = 0;

		while ( set->iter_ptr < set->size )
		{
			if ( !is_kmer_entity_null ( set->flags, set->iter_ptr ) )
			{
				node = set->array + set->iter_ptr;
				output1vt ( node, fp );
			}

			set->iter_ptr++;
		}
	}

	fprintf ( fp, "\n" );
	fprintf ( stderr, "%d vertex(es) output.\n", outvCounter );
	fclose ( fp );
	sprintf ( temp, "%s.preGraphBasic", outfile );
	fp = ckopen ( temp, "w" );
	fprintf ( fp, "VERTEX %d K %d\n", outvCounter, overlaplen );
	fprintf ( fp, "\nEDGEs %d\n", num_ed );
	fprintf ( fp, "\nMaxReadLen %d MinReadLen %d MaxNameLen %d\n", maxReadLen4all, minReadLen, maxNameLen );
	fclose ( fp );
}

void output_1edge ( preEDGE * edge, gzFile * fp )
{
	int i;
	gzprintf ( fp, ">length %d,", edge->length );
	print_kmer_gz ( fp, edge->from_node, ',' );
	print_kmer_gz ( fp, edge->to_node, ',' );
	gzprintf ( fp, "cvg %d, %d\n", edge->cvg, edge->bal_edge );

	for ( i = 0; i < edge->length; i++ )
	{
		gzprintf ( fp, "%c", int2base ( ( int ) edge->seq[i] ) );

		if ( ( i + 1 ) % 100 == 0 )
		{
			gzprintf ( fp, "\n" );
		}
	}

	if ( edge->length % 100 != 0 )
	{
		gzprintf ( fp, "\n" );
	}
}



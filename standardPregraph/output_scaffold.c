/*
 * output_scaffold.c
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

void output_contig_graph ( char * outfile )
{
	char name[256];
	FILE * fp;
	unsigned int i;
	sprintf ( name, "%s.contig.gvz", outfile );
	fp = ckopen ( name, "w" );
	fprintf ( fp, "digraph G{\n" );
	fprintf ( fp, "\tsize=\"512,512\";\n" );

	for ( i = num_ctg; i > 0; i-- )
	{
		fprintf ( fp, "\tV%d -> V%d[label =\"%d(%d)\"];\n", contig_array[i].from_vt, contig_array[i].to_vt, i, contig_array[i].length );
	}

	fprintf ( fp, "}\n" );
	fclose ( fp );
}
void output_scaf ( char * outfile )
{
	char name[256];
	FILE * fp;
	unsigned int i;
	CONNECT * connect;
	boolean flag;
	sprintf ( name, "%s.scaffold.gvz", outfile );
	fp = ckopen ( name, "w" );
	fprintf ( fp, "digraph G{\n" );
	fprintf ( fp, "\tsize=\"512,512\";\n" );

	for ( i = num_ctg; i > 0; i-- )
	{
		//if(contig_array[i].mask||!contig_array[i].downwardConnect)
		if ( !contig_array[i].downwardConnect )
		{
			continue;
		}

		connect = contig_array[i].downwardConnect;

		while ( connect )
		{
			//if(connect->mask||connect->deleted){
			if ( connect->deleted )
			{
				connect = connect->next;
				continue;
			}

			if ( connect->prevInScaf || connect->nextInScaf )
			{
				flag = 1;
			}
			else
			{
				flag = 0;
			}

			if ( !connect->mask )
				fprintf ( fp, "\tC%d_%d -> C%d_%d [label = \"%d(%d_%d)\"];\n", i, contig_array[i].length, connect->contigID, contig_array[connect->contigID].length,
				          connect->gapLen, flag, connect->weight );
			else
				fprintf ( fp, "\tC%d_%d -> C%d_%d [label = \"%d(%d_%d)\", color = red];\n", i, contig_array[i].length, connect->contigID, contig_array[connect->contigID].length,
				          connect->gapLen, flag, connect->weight );

			connect = connect->next;
		}
	}

	fprintf ( fp, "}\n" );
	fclose ( fp );
}

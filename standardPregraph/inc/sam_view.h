#ifndef SAM_VIEW_H
#define SAM_VIEW_H


static int g_min_mapQ = 0, g_flag_on = 0, g_flag_off = 0;
static char * g_library, *g_rg;
static int g_sol2sanger_tbl[128];

static void sol2sanger ( bam1_t * b )
{
	int l;
	uint8_t * qual = bam1_qual ( b );

	if ( g_sol2sanger_tbl[30] == 0 )
	{
		for ( l = 0; l != 128; ++l )
		{
			g_sol2sanger_tbl[l] = ( int ) ( 10.0 * log ( 1.0 + pow ( 10.0, ( l - 64 + 33 ) / 10.0 ) ) / log ( 10.0 ) + .499 );

			if ( g_sol2sanger_tbl[l] >= 93 ) { g_sol2sanger_tbl[l] = 93; }
		}
	}

	for ( l = 0; l < b->core.l_qseq; ++l )
	{
		int q = qual[l];

		if ( q > 127 ) { q = 127; }

		qual[l] = g_sol2sanger_tbl[q];
	}
}

static inline int __g_skip_aln ( const bam_header_t * h, const bam1_t * b )
{
	if ( b->core.qual < g_min_mapQ || ( ( b->core.flag & g_flag_on ) != g_flag_on ) || ( b->core.flag & g_flag_off ) )
		{ return 1; }

	if ( g_rg )
	{
		uint8_t * s = bam_aux_get ( b, "RG" );

		if ( s && strcmp ( g_rg, ( char * ) ( s + 1 ) ) == 0 ) { return 0; }
	}

	if ( g_library )
	{
		const char * p = bam_get_library ( ( bam_header_t * ) h, b );
		return ( p && strcmp ( p, g_library ) == 0 ) ? 0 : 1;
	}

	return 0;
}


#endif

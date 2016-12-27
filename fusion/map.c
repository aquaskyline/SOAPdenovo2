#include "stdinc.h" 
#include "newhash.h"
#include "extfunc.h" 
#include "extvab.h" 

//static void initenv(int argc, char **argv);


static void display_map_usage();

int call_align()
{
	time_t start_t,stop_t,time_bef,time_aft;
	time(&start_t);


	time(&time_bef);
	ctg_short = overlaplen+2;
	//printf("contig len cutoff: %d\n",ctg_short);
	prlContig2nodes(graphfile,ctg_short);
	time(&time_aft);
	//printf("time spent on De bruijn graph construction: %ds\n\n",
		//	(int)(time_aft-time_bef));
	//map read to edge one by one
	//printf("All contigs loaded");
	time(&time_bef);
	prlLongRead2Ctg(shortrdsfile,graphfile);
	time(&time_aft);
	//printf("time spent on mapping long reads: %ds\n\n",(int)(time_aft-time_bef));

	time(&time_bef);
	prlRead2Ctg(shortrdsfile,graphfile);
	time(&time_aft);
	//printf("time spent on mapping reads: %ds\n\n",(int)(time_aft-time_bef));
	
	free_Sets(KmerSets,thrd_num);

	time(&stop_t);
	//printf("overall time for alignment: %dm\n\n",(int)(stop_t-start_t)/60);
	printf("[%s]total time on mapping reads to contig :%dm\n",__FUNCTION__,(int)(stop_t-start_t)/60);
	return 0;
}

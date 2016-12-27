#include "stdinc.h" 
#include "newhash.h"
#include "extfunc.h" 
#include "extvab.h" 

static void initenv(int argc, char **argv);
static void display_scaff_usage();

static boolean LINK,SCAFF;


int call_scaffold()
{
	time_t start_t,stop_t,time_bef,time_aft;
	time(&start_t);

	//initenv(argc, argv);
	
	loadPEgrads(graphfile);
	
	time(&time_bef);
	loadUpdatedEdges(graphfile);
	time(&time_aft);
	//printf("time spent on loading edges %ds\n",(int)(time_aft-time_bef));

	if(!SCAFF){
		time(&time_bef);
		PE2Links(graphfile);
		time(&time_aft);
		//printf("time spent on loading pair end info %ds\n",(int)(time_aft-time_bef));
	
		time(&time_bef);
		Links2Scaf(graphfile);
		time(&time_aft);
		//printf("time spent on creating scaffolds %ds\n",(int)(time_aft-time_bef));

		scaffolding(100,graphfile);
	}

	prlReadsCloseGap(graphfile);


//	locateReadOnScaf(graphfile);

	free_pe_mem();
	if(index_array)
		free((void *)index_array);

	freeContig_array();
	
	//destroyPreArcMem();
	destroyConnectMem();
	deleteCntLookupTable();

	time(&stop_t);
	//printf("time elapsed: %dm\n",(int)(stop_t-start_t)/60);
	printf("[%s]total time on scaffolding : %d minute(s).\n",__FUNCTION__,(int)(stop_t-start_t)/60);

	return 0;
}

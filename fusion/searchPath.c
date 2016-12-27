#include "stdinc.h" 
#include "newhash.h"
#include "extfunc.h" 
#include "extvab.h" 

static int trace_limit = 5000; //the times function is called in a search
/*
 search connection paths which were masked along related contigs
 start from one contig, end with another
 path length includes the length of the last contig
*/
void traceAlongMaskedCnt(unsigned int destE,unsigned int currE,int max_steps,int min,int max,
		int index,int len,int *num_route)
{
	num_trace++;
	if(num_trace>trace_limit||*num_route>=max_n_routes){
		return;
	}

	unsigned int *array;
	int num,i,length;
	CONNECT *ite_cnt;

	if(index>0)// there're at most max_steps edges stored in this array including the destination edge
		length = len + contig_array[currE].length; 
	else 
		length = 0;
	if(index>max_steps||length>max)
		return;  // this is the only situation we stop 
	if(index>0)// there're at most max_steps edges stored in this array including the destination edge
		so_far[index-1] = currE;  

	if(currE==destE&&index==0){
		printf("traceAlongMaskedCnt: start and destination are the same\n");
		return;
	}

	if(currE==destE && length>=min &&length<=max){
		num = *num_route;
		array = found_routes[num];
		for(i=0;i<index;i++)
			array[i] = so_far[i];
		if(index<max_steps)   
			array[index] = 0;  //indicate the end of the route
		
		*num_route = ++num;
	}   // one route is extrated, but we don't terminate searching

	ite_cnt = contig_array[currE].downwardConnect;
	while(ite_cnt){
		if(!ite_cnt->mask||ite_cnt->deleted){
			ite_cnt = ite_cnt->next;
			continue;
		}
		traceAlongMaskedCnt(destE,ite_cnt->contigID,max_steps,min,max,
			index+1,length + ite_cnt->gapLen,num_route);
		ite_cnt = ite_cnt->next;	
	}

}
// search connection paths from one connect to a contig
// path length includes the length of the last contig
void traceAlongConnect(unsigned int destE,CONNECT *currCNT,int max_steps,int min,int max,int index,int len,int *num_route)
{
	num_trace++;
	if(num_trace>trace_limit||*num_route>=max_n_routes){
		return;
	}

	unsigned int *array,currE;
	int num,i,length;
	CONNECT *ite_cnt;

	currE = currCNT->contigID;
	length = len + currCNT->gapLen;
	length += contig_array[currE].length; 
	
	if(index>max_steps||length>max)
		return;  // this is the only situation we stop 
		/*
	if(globalFlag)
		printf("B: step %d, ctg %d, length %d\n",index,currCNT->contigID,length);
		*/
	if(currE==destE&&index==1){
		printf("traceAlongConnect: start and destination are the same\n");
		return;
	}
	
	so_far[index-1] = currE;  // there're at most max_steps edges stored in this array including the destination edge

	if(currE==destE && length>=min &&length<=max){
		num = *num_route;
		array = found_routes[num];
		for(i=0;i<index;i++)
			array[i] = so_far[i];
		if(index<max_steps)   
			array[index] = 0;  //indicate the end of the route
		
		*num_route = ++num;
	}   // one route is extrated, but we don't terminate searching

	if(currCNT->nextInScaf){
		traceAlongConnect(destE,currCNT->nextInScaf,max_steps,min,max,index+1,length,num_route);
		return;
	}

	ite_cnt = contig_array[currE].downwardConnect;
	while(ite_cnt){
		if(ite_cnt->mask||ite_cnt->deleted){
			ite_cnt = ite_cnt->next;
			continue;
		}
		traceAlongConnect(destE,ite_cnt,max_steps,min,max,index+1,length,num_route);
		ite_cnt = ite_cnt->next;	
	}

}

//find paths in the graph from currE to destE, its length does not include length of both end contigs
void traceAlongArc(unsigned int destE,unsigned int currE,int max_steps,int min,int max,int index,int len,int *num_route)
{
	num_trace++;
	if(num_trace>trace_limit||*num_route>=max_n_routes){
		return;
	}

	unsigned int *array,out_ed,vt;
	int num,i,pos,length;
	preARC *parc;

	pos = index;
	if(pos>max_steps||len>max)
		return;  // this is the only situation we stop 
	if(currE==destE&&pos==0){
		printf("traceAlongArc: start and destination are the same\n");
		return;
	}
	
	if(pos>0)  // pos starts with 0 for the starting edge
		so_far[pos-1] = currE;         // there're at most max_steps edges stored in this array including the destination edge

	if(currE==destE && len>=min){
		num = *num_route;
		array = found_routes[num];
		for(i=0;i<pos;i++)
			array[i] = so_far[i];
		if(pos<max_steps)   
			array[pos] = 0;  //indicate the end of the route
		
		*num_route = ++num;
	}   // one route is extrated, but we don't terminate searching
	if(pos==max_steps||len==max)
		return;
	if(pos++>0)  //not the starting edge
		length = len + contig_array[currE].length;
	else
		length = len;

	
	vt = contig_array[currE].to_vt;
	
	parc = contig_array[currE].arcs;
	while(parc){
		out_ed = parc->to_ed;
		traceAlongArc(destE,out_ed,max_steps,min,max,pos,length,num_route);
		parc = parc->next;	
	}

}

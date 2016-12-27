#include "stdinc.h" 
#include "newhash.h"
#include "extfunc.h" 
#include "extvab.h" 
#include "dfibHeap.h"
#include "fibHeap.h"
#include "darray.h"


//static CTGinHEAP *ctg4heapArray;
extern int inputLinks(FILE *fp, int insertS,char *line);
//unsigned int traverse(unsigned int node,int *far_count,unsigned int *farpath,
	//int *curr_count,unsigned int *currpath,int *used_count,unsigned int *used,int *max_dist,int *node_dist);
//static int *sub_arr;
//static int sub_counter=0;
int rev_comp (const void * a, const void * b)
{
  return ( *(int*)b - *(int*)a );
}
void potential()
{

	char name[256],*line;
	FILE *fp;
	int i;
	int flag2;

	loadUpdatedEdges(graphfile);

	if(!pes)
		loadPEgrads(graphfile);
	
	sprintf(name,"%s.links",graphfile);
	fp = ckopen(name,"r");
	
	createCntMemManager();
	createCntLookupTable();

	lineLen = 1024;
	line = (char *)ckalloc(lineLen*sizeof(char));

	fgets(line,lineLen,fp);
	line[0] = '\0';
	fprintf(stderr,"[%s]before inputLinks loop.\n",__FUNCTION__);
	for(i=0;i<gradsCounter;i++){	
		flag2 = inputLinks(fp,pes[i].insertS,line);
	}
	fprintf(stderr,"[%s]links file loaded.\n",__FUNCTION__);
	//ctg4heapArray=ckalloc(100000*sizeof(CTGinHEAP));
	//unsigned int *farthest_path=ckalloc(1000000*sizeof(unsigned int));
	//int farthest_boarder;
	//STACK *track=createStack(100000,sizeof(unsigned int));
	unsigned int *curr_path=ckalloc(1000000*sizeof(unsigned int));
	int curr_boarder;
	unsigned int *dist=ckalloc(10000000*sizeof(unsigned int));
	//int dist_boarder;
	int *predict=ckalloc(1000000*sizeof(int));
	int pred_count=0;
	int used=0;
	//CONNECT *cnt_stack=ckalloc(10000000*sizeof(CONNECT*));
	//int cnt_count=0;
	for(i=1;i<=num_ctg;i++){
		if(contig_array[i].inSubGraph)
			continue;
		if(!contig_array[i].downwardConnect){
			predict[pred_count++]=contig_array[i].length;
			contig_array[i].inSubGraph=1;
			contig_array[getTwinCtg(i)].inSubGraph=1;
			fprintf(stderr,"[%d] traversed %d %d .\n",__LINE__,i,getTwinCtg(i));
			++used;
			continue;
		}
		//depth first traversal
		//farthest_boarder=0;
		curr_boarder=0;
		//used_boarder=0;
		int max_dist=0;
		//int node_dist=0;
		int len=0;
		//contig_array[i].inSubGraph=1;
		//contig_array[getTwinCtg(i)].inSubGraph=1;
		//if(contig_array[i].downwardConnect){
		//traverse(i,&farthest_boarder,farthest_path,&curr_boarder,curr_path,&used_boarder,used,&max_dist,&node_dist);
		//cnt_stack[0]=contig_array[i].downwardConnect;//put in stack
		//++cnt_count;
		
		contig_array[i].inSubGraph=1;
		contig_array[getTwinCtg(i)].inSubGraph=1;
		fprintf(stderr,"[%d] traversed %d %d .\n",__LINE__,i,getTwinCtg(i));
		++used;
		curr_path[curr_boarder]=i;
		dist[curr_boarder]=0;
		while(curr_boarder>=0){
			int curr_bound=curr_boarder;
			int curr_node=curr_path[curr_boarder--];
			int base_dist=dist[curr_bound];
			CONNECT *curr_cnt=contig_array[curr_node].downwardConnect;
			while(curr_cnt){//push all adjacent connect
				if(curr_cnt->weight<3||contig_array[curr_cnt->contigID].inSubGraph
					||contig_array[getTwinCtg(curr_cnt->contigID)].inSubGraph){
					curr_cnt=curr_cnt->next;
					continue;
				}
				curr_path[++curr_boarder]=curr_cnt->contigID;
				contig_array[curr_cnt->contigID].inSubGraph=1;
				contig_array[getTwinCtg(curr_cnt->contigID)].inSubGraph=1;
				fprintf(stderr,"[%d] traversed %d %d .\n",__LINE__,curr_cnt->contigID,getTwinCtg(curr_cnt->contigID));
				++used;
				dist[curr_boarder]=base_dist+
					curr_cnt->gapLen+contig_array[curr_cnt->contigID].length;
				
				if(dist[curr_boarder]>max_dist)
					max_dist=dist[curr_boarder];
				//fprintf(stderr,"curr_boarder %d node_dist %d max_dist %d \n",curr_boarder,
				//	dist[curr_boarder],max_dist);
				curr_cnt=curr_cnt->next;
			}
			
		}
		len+=max_dist;
		
		//}
		if(contig_array[getTwinCtg(i)].downwardConnect){
			curr_boarder=0;
			curr_path[curr_boarder]=i;
			dist[curr_boarder]=0;
			
			while(curr_boarder>=0){
				int curr_bound=curr_boarder;
				int curr_node=curr_path[curr_boarder--];
				int base_dist=dist[curr_bound];
				CONNECT *curr_cnt=contig_array[curr_node].downwardConnect;
				while(curr_cnt){//push all adjacent connect
					if(curr_cnt->weight<3||contig_array[curr_cnt->contigID].inSubGraph
						||contig_array[getTwinCtg(curr_cnt->contigID)].inSubGraph){
						curr_cnt=curr_cnt->next;
						continue;
					}
					curr_path[++curr_boarder]=curr_cnt->contigID;
					contig_array[curr_cnt->contigID].inSubGraph=1;
					contig_array[getTwinCtg(curr_cnt->contigID)].inSubGraph=1;
					fprintf(stderr,"[%d] traversed %d %d .\n",__LINE__,curr_cnt->contigID,getTwinCtg(curr_cnt->contigID));
					++used;
					dist[curr_boarder]=base_dist+
						curr_cnt->gapLen+contig_array[curr_cnt->contigID].length;
					
					if(dist[curr_boarder]>max_dist)
						max_dist=dist[curr_boarder];
					//fprintf(stderr,"curr_boarder %d node_dist %d max_dist %d \n",curr_boarder,
					//	dist[curr_boarder],max_dist);
					curr_cnt=curr_cnt->next;
				}
				
			}
			len+=max_dist;
		}
		/*int ii;
		for(ii=0;ii<used_boarder;++ii){
			contig_array[used[ii]].inSubGraph=0;
		}*/
		if(len!=0){
			predict[pred_count++]=len;
			fprintf(stderr,"[%s]contig %d effective with length %d.\n",__FUNCTION__,i,len);
		}
		fprintf(stderr,"[%s]contig %d over.\n",__FUNCTION__,i);
		
		
	}

	free((void *)line);
	fclose(fp);
	long long int sum=0;
	for(i=0;i<pred_count;++i){
		sum+=predict[i];
	}
	long long int half=sum/2;
	printf("sum %lld , half  %lld.\n",sum,half);
	qsort(predict,pred_count,sizeof(int),rev_comp);
	sum=0;
	for(i=0;i<pred_count;++i){	
		printf("len:\t%d\n",predict[i]);
	}
	for(i=0;i<pred_count;++i){
		sum+=predict[i];
		printf("len:\t%d\n",predict[i]);
		if(sum>=half)
			break;
	}
	printf("N50 %d , half %lld.\n",predict[i],half);
	printf("used contig %d",used);
}

/*
unsigned int traverse(unsigned int node,int *far_count,unsigned int *farpath,
	int *curr_count,unsigned int *currpath,int *used_count,unsigned int *used,int *max_dist,int *node_dist){
	unsigned int bal = getTwinCtg(node);
	
	currpath[(*curr_count)++]=node;
	used[(*used_count)++]=node;
	used[(*used_count)++]=bal;
	contig_array[node].inSubGraph=1;
	contig_array[bal].inSubGraph=1;
	
	fprintf(stderr,"farcount %d curr_count %d node_dist %d max_dist %d.\n",*far_count,*curr_count,*node_dist,*max_dist);
	CONNECT *tmp_cnt=contig_array[node].downwardConnect;
	while(tmp_cnt){
		unsigned int ctg,bal_ctg;
		ctg=tmp_cnt->contigID;
		bal_ctg=getTwinCtg(ctg);
		if(contig_array[ctg].inSubGraph||contig_array[bal_ctg].inSubGraph
			||contig_array[ctg].flag||contig_array[bal_ctg].flag){
			tmp_cnt=tmp_cnt->next;
			continue;	
		}
		*node_dist+=(tmp_cnt->gapLen+contig_array[ctg].length);
		if(*node_dist>*max_dist){
			int i;
			for(i=0;i<*curr_count;++i){
				farpath[i]=currpath[i];
			}
			*far_count=*curr_count;
			*max_dist=*node_dist+tmp_cnt->gapLen;
		}
		traverse(tmp_cnt->contigID,far_count,farpath,curr_count,currpath,used_count,used,max_dist,node_dist);
		*node_dist-=(tmp_cnt->gapLen+contig_array[ctg].length);
		tmp_cnt=tmp_cnt->next;
	}
	--(*curr_count);
	
	return 0;
}
*/

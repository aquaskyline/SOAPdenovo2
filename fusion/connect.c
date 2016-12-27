#include "stdinc.h" 
#include "newhash.h"
#include "extfunc.h" 
#include "extvab.h" 

#define CNBLOCKSIZE 100000

void createCntMemManager()
{
	if(!cn_mem_manager)
		cn_mem_manager = createMem_manager(CNBLOCKSIZE,sizeof(CONNECT));
	//else
		//printf("cn_mem_manger was created\n");
}

void destroyConnectMem()
{
	freeMem_manager(cn_mem_manager);
	cn_mem_manager = NULL;
}

CONNECT *allocateCN(unsigned int contigId, int gap)
{
	CONNECT *newCN;
	newCN = (CONNECT *)getItem(cn_mem_manager);
	newCN->contigID = contigId;
	newCN->gapLen = gap;
	
	newCN->minGap = 0;
	newCN->maxGap = 0;
	newCN->bySmall = 0;
	newCN->weakPoint = 0;

	newCN->weight = 1;
	newCN->weightNotInherit = 0;
	newCN->mask = 0;
	newCN->used = 0;
	newCN->checking = 0;
	newCN->deleted = 0;
	newCN->prevInScaf = 0;
	newCN->inherit = 0;
	newCN->singleInScaf = 0;
	newCN->nextInScaf = NULL;
	newCN->PE=NULL;//(int *)ckalloc(CNBLOCKSIZE*sizeof(int));
	
	return newCN;
}

void output_cntGVZ(char *outfile)
{	
	char name[256];
	FILE *fp;
	unsigned int i;
	CONNECT *connect;
	boolean flag;

	sprintf(name,"%s.scaffold.gvz",outfile);
	fp = ckopen(name,"w");
	fprintf(fp,"digraph G{\n");
	fprintf(fp,"\tsize=\"512,512\";\n");

	for(i=num_ctg;i>0;i--){
		//if(contig_array[i].mask||!contig_array[i].downwardConnect)
		if(!contig_array[i].downwardConnect)
			continue;
		connect = contig_array[i].downwardConnect;
		while(connect){
			//if(connect->mask||connect->deleted){
			if(connect->deleted){
				connect = connect->next;
				continue;
			}
			if(connect->prevInScaf||connect->nextInScaf)
				flag = 1;
			else
				flag = 0;
			if(!connect->mask)
				fprintf(fp,"\tC%d_%d -> C%d_%d [label = \"%d(%d_%d)\"];\n"
				,i,contig_array[i].length,connect->contigID,contig_array[connect->contigID].length,
				connect->gapLen,flag,connect->weight);
			else	
				fprintf(fp,"\tC%d_%d -> C%d_%d [label = \"%d(%d_%d)\", color = red];\n"
				,i,contig_array[i].length,connect->contigID,contig_array[connect->contigID].length,
				connect->gapLen,flag,connect->weight);
			connect = connect->next;
		}
	}
	fprintf(fp,"}\n");
	fclose(fp);
}

/***************** below this line all codes are about lookup table *****************/

void createCntLookupTable()
{
	if(!cntLookupTable)
		cntLookupTable = (CONNECT **)ckalloc((3*num_ctg+1)*sizeof(CONNECT *));
}

void deleteCntLookupTable()
{
	if(cntLookupTable){
		free((void *)cntLookupTable);
		cntLookupTable = NULL;
	}
}

void putCnt2LookupTable(unsigned int from_c,CONNECT *cnt)
{
	if(!cnt||!cntLookupTable)
		return;
	unsigned int index = 2*from_c + cnt->contigID;
	cnt->nextInLookupTable = cntLookupTable[index];
	cntLookupTable[index] = cnt;
}

static CONNECT *getCntInLookupTable(unsigned int from_c,unsigned int to_c)
{
	unsigned int index = 2*from_c + to_c;
	CONNECT *ite_cnt = cntLookupTable[index];
	while(ite_cnt){
		if(ite_cnt->contigID==to_c)
			return ite_cnt;
		ite_cnt = ite_cnt->nextInLookupTable;	
	}
	return NULL;
}

CONNECT *getCntBetween(unsigned int from_c, unsigned int to_c)
{
	CONNECT *pcnt;
		
	if(cntLookupTable){
		pcnt = getCntInLookupTable(from_c,to_c);
		return pcnt;
	}
	pcnt = contig_array[from_c].downwardConnect;

	while(pcnt){
		if(pcnt->contigID==to_c)
			return pcnt;
		pcnt = pcnt->next;	
	}
	return pcnt;
}
/*
void removeCntInLookupTable(unsigned int from_c,unsigned int to_c)
{
	unsigned int index = 2*from_c + to_c;
	CONNECT *ite_cnt = cntLookupTable[index];
	CONNECT *cnt;

	if(!ite_cnt){
		printf("removeCntInLookupTable: not found A\n");
		return;
	}
	if(ite_cnt->contigID==to_c){
		cntLookupTable[index] = ite_cnt->nextInLookupTable;	
		return;
	}
	
	while(ite_cnt->nextInLookupTable&&ite_cnt->nextInLookupTable->contigID!=to_c)
		ite_cnt = ite_cnt->nextInLookupTable;	
	
	if(ite_cnt->nextInLookupTable){
		cnt = ite_cnt->nextInLookupTable;
		ite_cnt->nextInLookupTable = cnt->nextInLookupTable;
		return;
	}
	printf("removeCntInLookupTable: not found B\n");
	return;
}
*/

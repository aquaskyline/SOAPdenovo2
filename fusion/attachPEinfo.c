#include "stdinc.h" 
#include "newhash.h"
#include "extfunc.h" 
#include "extvab.h" 
#include "stack.h"

#define CNBLOCKSIZE 10000
#define GAPARRSIZE 256
#define BIG_NEG -10000000
#define BIG_POS 10000000
static STACK * isStack;
static int ignorePE1,ignorePE2,ignorePE3,ignorePE4,ignorePE5,static_flag;
static int onsameCtgPE;
static unsigned long long peSUM;

//static boolean staticF;

static int existCounter;

int calcuIS(STACK *intStack,int *SD);


static int cmp_pe(const void *a,const void *b)
{
	PE_INFO *A,*B;
	A = (PE_INFO *)a;
	B = (PE_INFO *)b;

	if(A->rank>B->rank)
		return 1;
	else if(A->rank==B->rank)
		return 0;
	else
		return -1;
}

void loadPEgrads(char *infile)
{
	FILE *fp;
	char name[256],line[1024];
	int i;
	boolean rankSet=1;

	sprintf(name,"%s.peGrads",infile);
    	fp = fopen(name,"r");
	if(!fp){
		printf("can not open file %s .\n",name);
		gradsCounter = 0;
		return;
	}
	
	while(fgets(line,sizeof(line),fp)!=NULL){
		if(line[0] == 'g'){
			sscanf(line+10, "%d %lld %d",&gradsCounter,&n_solexa,&maxReadLen);
			//printf("there're %d grads, %lld reads, max read len %d\n",gradsCounter,n_solexa,maxReadLen);
			printf("[%s]reads statistic : %lld reads with max len %d in %d grads .\n",__FUNCTION__,n_solexa,maxReadLen,gradsCounter);
			break;
		}
	}

	alloc_pe_mem(gradsCounter);

	for(i=0;i<gradsCounter;i++){
		fgets(line,sizeof(line),fp);
		pes[i].rank = 0;
		sscanf(line,"%d %lld %d %d",&(pes[i].insertS),&(pes[i].PE_bound),&(pes[i].rank),&(pes[i].pair_num_cut));
		if(pes[i].rank<1)
			rankSet = 0;
	}
	fclose(fp);
	if(rankSet){
		qsort(&pes[0],gradsCounter,sizeof(PE_INFO),cmp_pe);
		return;
	}
	int lastRank = 0;
	for(i=0;i<gradsCounter;i++){
		if(i==0)
			pes[i].rank = ++lastRank;
		else if(pes[i].insertS<300)
			pes[i].rank = lastRank;
		else if(pes[i].insertS<800){
			if(pes[i-1].insertS<300)
				pes[i].rank = ++lastRank;
			else
				pes[i].rank = lastRank;
		}else if(pes[i].insertS<3000){
			if(pes[i-1].insertS<800)
				pes[i].rank = ++lastRank;
			else
				pes[i].rank = lastRank;
		}else if(pes[i].insertS<7000){
			if(pes[i-1].insertS<3000)
				pes[i].rank = ++lastRank;
			else
				pes[i].rank = lastRank;
		}else{
			if(pes[i-1].insertS<7000)
				pes[i].rank = ++lastRank;
			else
				pes[i].rank = lastRank;
		}
	}
		
}


CONNECT *add1Connect(unsigned int e1, unsigned int e2, int gap, int weight,boolean inherit)
{
	if(e1==e2||e1==getTwinCtg(e2))
		return NULL;
	CONNECT *connect=NULL;
	long long sum;
	if(weight>255)
		weight = 255;

	connect = getCntBetween(e1, e2);
	if(connect){
		if(!weight)
			return connect;
		existCounter++;
		if(!inherit){
			sum = connect->weightNotInherit*connect->gapLen + gap*weight;
			connect->gapLen = sum/(connect->weightNotInherit+weight);
			if(connect->weightNotInherit+weight <=255)
				connect->weightNotInherit += weight;
			else if(connect->weightNotInherit<255)
				connect->weightNotInherit = 255;
		}else{
			sum = connect->weight*connect->gapLen + gap*weight;
			connect->gapLen = sum/(connect->weight+weight);
			if(!connect->inherit){
				connect->maxSingleWeight = connect->weightNotInherit;
			}
			connect->inherit = 1;
			connect->maxSingleWeight = connect->maxSingleWeight>weight ?
					connect->maxSingleWeight:weight;
		}
		if(connect->weight+weight <=255){
			connect->weight += weight;
		}else if(connect->weight<255){
			connect->weight = 255;
		}
		
	}else{
		newCntCounter++;
		connect = allocateCN(e2,gap);
		if(cntLookupTable)
			putCnt2LookupTable(e1,connect);
		connect->weight = weight;
		if(contig_array[e1].mask||contig_array[e2].mask){
			connect->mask = 1;
		}
		connect->next = contig_array[e1].downwardConnect;
		contig_array[e1].downwardConnect = connect;
		if(!inherit){
			connect->weightNotInherit = weight;
		}else{
			connect->weightNotInherit = 0;
			connect->inherit = 1;
			connect->maxSingleWeight = weight;
		}	
	}
	
	return connect;
}
CONNECT *add1AccuConnect(unsigned int e1, unsigned int e2, int gap, int weight)
{
	if(e1==e2||e1==getTwinCtg(e2))
		return NULL;
	CONNECT *connect=NULL;
	//long long sum;
	if(weight>255)
		weight = 255;

	connect = getCntBetween(e1, e2);
	if(connect){
		if(!weight)
			return connect;
		existCounter++;
		//if(!inherit){
			//sum = connect->weightNotInherit*connect->gapLen + gap*weight;
			//connect->gapLen = sum/(connect->weightNotInherit+weight);
			int i=connect->weightNotInherit;
			
			if(connect->weightNotInherit+weight <=255)
				connect->weightNotInherit += weight;
			else if(connect->weightNotInherit<255)
				connect->weightNotInherit = 255;
			for(;i<connect->weightNotInherit;i++){
				connect->PE[i]=gap;
				fprintf(stderr,"inputting a PE with estimated gap size %d\n",gap);
			}
		/*}else{
			//sum = connect->weight*connect->gapLen + gap*weight;
			//connect->gapLen = sum/(connect->weight+weight);
			if(!connect->inherit){
				connect->maxSingleWeight = connect->weightNotInherit;
			}
			connect->inherit = 1;
			connect->maxSingleWeight = connect->maxSingleWeight>weight ?
					connect->maxSingleWeight:weight;
		}*/
		if(connect->weight+weight <=255){
			connect->weight += weight;
		}else if(connect->weight<255){
			connect->weight = 255;
		}
		
	}else{
		newCntCounter++;
		connect = allocateCN(e2,gap);
		if(cntLookupTable)
			putCnt2LookupTable(e1,connect);
		connect->weight = weight;
		connect->PE=(int *)ckalloc(GAPARRSIZE*sizeof(int));//newly added
		fprintf(stderr,"creating array for PEs in a connection.\n");
		int i;
		for(i=0;i<weight;i++){
			connect->PE[i]=gap;
			fprintf(stderr,"inputting a PE with estimated gap size %d\n",gap);
		}
		if(contig_array[e1].mask||contig_array[e2].mask){
			connect->mask = 1;
		}
		connect->next = contig_array[e1].downwardConnect;
		contig_array[e1].downwardConnect = connect;
		//if(!inherit){
			connect->weightNotInherit = weight;
		/*}else{
			connect->weightNotInherit = 0;
			connect->inherit = 1;
			connect->maxSingleWeight = weight;
		}*/	
	}
	
	return connect;
}
int attach1PE(unsigned int e1,int pre_pos,unsigned int bal_e2,int pos,int insert_size)
{
	int gap,realpeSize;
	unsigned int bal_e1,e2;
	if(e1==bal_e2){
		ignorePE1++;
		return -1;  //orientation wrong
	}

	bal_e1 = getTwinCtg(e1);
	e2 = getTwinCtg(bal_e2);
	if(e1==e2){
		realpeSize = contig_array[e1].length + overlaplen - pre_pos - pos;
		if(realpeSize>0){
			peSUM += realpeSize;
			onsameCtgPE++;
			if((int)contig_array[e1].length>insert_size){
				int *item = (int *)stackPush(isStack);
				(*item) = realpeSize;
			}
		}
		return 2;
	}

	gap = insert_size - overlaplen + pre_pos + pos - contig_array[e1].length - contig_array[e2].length;
	//fprintf(stderr,"[%s]\tgap\t%d\t%d\t%f\t%f\n",__FUNCTION__,gap,insert_size,close_threshold,insert_size*close_threshold);
	if(gap<-(insert_size*close_threshold)){
		ignorePE2++;
		return 0;
	}
	if(gap>insert_size){
		ignorePE3++;
		return 0;
	}
	add1AccuConnect(e1,e2,gap,1);
	add1AccuConnect(bal_e2,bal_e1,gap,1);

	return 1;
}

int connectByPE_grad(FILE *fp,int peGrad,char *line)
{
	fprintf(stderr,"[%s]entering this function.\n",__FUNCTION__);
	long long pre_readno,readno,minno,maxno;
	int pre_pos,pos,flag,PE,count=0;
	unsigned int pre_contigno,contigno,newIndex;

	if(peGrad<0||peGrad>gradsCounter){
		printf("[%s]specified pe grad is out of bound .\n",__FUNCTION__);
		return 0;
	}
	maxno = pes[peGrad].PE_bound;
	if(peGrad==0)
		minno = 0;
	else
		minno = pes[peGrad-1].PE_bound;
		
	onsameCtgPE = peSUM = 0;
	PE = pes[peGrad].insertS;
	if(strlen(line)){
			sscanf(line,"%lld %d %d",&pre_readno,&pre_contigno,&pre_pos);
			//printf("first record %d %d %d\n",pre_readno,pre_contigno,pre_pos);
			if(pre_readno<=minno)
				pre_readno = -1;
	}
	else
		pre_readno = -1;
	ignorePE1 = ignorePE2 = ignorePE3 = ignorePE4 = ignorePE5 = 0;
	static_flag = 1;
	isStack = (STACK *)createStack(CNBLOCKSIZE,sizeof(int));
	while(fgets(line,lineLen,fp)!=NULL){
			sscanf(line,"%lld %d %d",&readno,&contigno,&pos);
			if(readno>maxno)
				break;
			if(readno<=minno)
				continue;
			
			newIndex = index_array[contigno];
			//if(contig_array[newIndex].bal_edge==0)
			if(isSameAsTwin(newIndex))
				continue;
			if(PE&&(readno%2==0)&&(pre_readno==readno-1)){ // they are a pair of reads
				flag = attach1PE(pre_contigno,pre_pos,newIndex,pos,PE);
				if(flag==1)
					count++;
			}
			pre_readno = readno;
			pre_contigno = newIndex;
			pre_pos = pos;
	}
	printf("[%s]Finish loading all PEs in grad %d .\n",__FUNCTION__,peGrad);
	printf("[%s]Calculating estimated gap size for all connections .\n",__FUNCTION__);
	unsigned int i;
	for(i=1;i<=num_ctg;i++){
			CONNECT *tmp=contig_array[i].downwardConnect;
			while(tmp){
				if(tmp->weightNotInherit<=8&&tmp->weightNotInherit>2){//delete max and min value
					int max=BIG_NEG,maxid=-1,min=BIG_POS,minid=-1;
					int weight=tmp->weightNotInherit;
					int ii;
					for(ii=0;ii<weight;ii++){
						if(tmp->PE[ii]>max){
							max=tmp->PE[ii];
							maxid=ii;
						}
						if(tmp->PE[ii]<=min){
							min=tmp->PE[ii];
							minid=ii;
						}
					}
					int sum=0;
					for(ii=0;ii<weight;ii++){
						if(ii!=maxid&&ii!=minid){
							sum+=tmp->PE[ii];
						}
					}
					ignorePE4+=2;
					tmp->gapLen=sum/(weight-2);
					fprintf(stderr,"estimating contigs' gap by removing max&min PE ,with max&min %d %d\n",
						tmp->PE[maxid],tmp->PE[minid]);
				}else if(tmp->weightNotInherit>8){//delete values exceed 3*SD
					long long int sum=0;
					int weight=tmp->weightNotInherit;
					int ii;
					int counter=0;
					for(ii=0;ii<weight;ii++){
						sum+=tmp->PE[ii];
					}
					
					long long int avg=sum/weight;
					sum = 0;
					for(ii=0;ii<weight;ii++){
						sum+=((tmp->PE[ii]-avg)*(tmp->PE[ii]-avg));
					}
					
					double SD=(sqrt((double)sum/(weight-1)))*3;//just for fast
					sum=0;
					int num=0;
					for(ii=0;ii<weight;ii++){
						if(abs(tmp->PE[ii]-avg)<=SD){
							sum+=tmp->PE[ii];
							num++;
						}else{
							ignorePE5++;
							counter++;
						}
					}
					if(num==0){
						fprintf(stderr,"[%s]num=0 in removing exceed 3*SD(%.1f) avg(%d)step",__FUNCTION__,SD,avg);
						for(ii=0;ii<weight;ii++){
							fprintf(stderr,"%d\t",tmp->PE[ii]);
						}
					}
					tmp->gapLen=sum/num;
					fprintf(stderr,"estimating contigs' gap by removing PE exceeding 3*SD ,removing %d PEs\n",counter);
				}else if(tmp->weightNotInherit<=2){
					int weight=tmp->weightNotInherit;
					int sum=0;
					int ii;
					for(ii=0;ii<weight;ii++){
						sum+=tmp->PE[ii];
					}
					tmp->gapLen=sum/weight;
					fprintf(stderr,"weight too small , directly estimate gap size.\n");
				}
				//fprintf(stderr,"finish %d connection.\n",i);
				free((void *)tmp->PE);
				tmp=tmp->next;
			}
	}
	//printf("%d PEs with insert size %d attached, %d + %d + %d ignored\n",count,PE,ignorePE1,ignorePE2,ignorePE3);
	fprintf(stderr,"[%s]%d PEs of insert size %d loaded .\n",__FUNCTION__,count,PE);
	fprintf(stderr,"[%s]PEs discarded:%d because of wrong orientation,%d too close,%d too far,\n",__FUNCTION__,ignorePE1,ignorePE2,ignorePE3);
	fprintf(stderr,"[%s]%d deleted by removing max&min , %d not fall in 3*SD.\n",__FUNCTION__,ignorePE4,ignorePE5);
	printf("[%s]%d PEs of insert size %d loaded .\n",__FUNCTION__,count,PE);
	printf("[%s]PEs discarded :%d because of wrong orientation,%d too close,%d too far ,\n",__FUNCTION__,ignorePE1,ignorePE2,ignorePE3);
	printf("[%s]%d deleted by removing max&min , %d not fall in 3*SD .\n",__FUNCTION__,ignorePE4,ignorePE5);
	
	if(onsameCtgPE>0){
		//printf("estimated PE size %lli, by %d pairs\n",peSUM/onsameCtgPE,onsameCtgPE);
		int SD=0;
		int avg=calcuIS(isStack,&SD);
		printf("[%s]%d PE attached on same contig with estimated insert size %d SD %d .\n",__FUNCTION__,onsameCtgPE,avg,SD);
	}
	//printf("on contigs longer than %d, %d pairs found,",PE,isStack->item_c);
	//printf("insert_size estimated: %d\n",calcuIS(isStack));
	freeStack(isStack);
	return count;
}


int calcuIS(STACK *intStack,int *SD)
{
	long long sum=0;
	int avg=0;
	int *item;
	int num = intStack->item_c;
	
	if(num<100)
		return avg;
	stackBackup(intStack);
	while((item=(int *)stackPop(intStack))!=NULL)
		sum += *item;

	stackRecover(intStack);
	num = intStack->item_c;
	avg = sum/num;

	sum = 0;
	stackBackup(intStack);
	while((item=(int *)stackPop(intStack))!=NULL)
		sum += (*item-avg)*(*item-avg);
	
	*SD = sqrt(sum/(num-1));
	if(SD==0){
		//printf("SD=%d, ",SD);
		return avg;
	}
	stackRecover(intStack);
	sum = num = 0;
	while((item=(int *)stackPop(intStack))!=NULL)
		if(abs(*item-avg)<3**SD){
			sum += *item;
			num++;
		}
	
	avg = sum/num;
	//printf("SD=%d, ",SD);
	return avg;

}

unsigned int getTwinCtg(unsigned int ctg)
{
	return ctg + contig_array[ctg].bal_edge - 1;
}

boolean isSmallerThanTwin(unsigned int ctg)
{
	return contig_array[ctg].bal_edge > 1;
}

boolean isLargerThanTwin(unsigned int ctg)
{
	return contig_array[ctg].bal_edge < 1;
}

boolean isSameAsTwin(unsigned int ctg)
{
	return contig_array[ctg].bal_edge == 1;
}

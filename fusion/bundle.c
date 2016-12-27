#include "stdinc.h" 
#include "newhash.h"
#include "extfunc.h" 
#include "extvab.h" 
#include "dfibHeap.h"
#include "fibHeap.h"
#include "darray.h"


#define CNBLOCKSIZE 10000
#define GAPARRSIZE 256
#define BIG_NEG -10000000
#define BIG_POS 10000000
static STACK * isStack;
static int onsameCtgPE;
extern int calcuIS(STACK *intStack,int *SD);
void outputBundle(FILE *fp, int insertS);

static CONNECT *bun1AccuConnect(unsigned int e1, unsigned int e2, int gap, int weight)
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
		//existCounter++;
		//if(!inherit){
			//sum = connect->weightNotInherit*connect->gapLen + gap*weight;
			//connect->gapLen = sum/(connect->weightNotInherit+weight);
			int i=connect->weightNotInherit;
			
			if(connect->weightNotInherit+weight <=255)
				connect->weightNotInherit += weight;
			else if(connect->weightNotInherit<255)
				connect->weightNotInherit = 255;
			for(;i<connect->weightNotInherit;i++){
				//connect->PE[i]=gap;
				//fprintf(stderr,"inputting a PE with estimated gap size %d\n",gap);
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
		//connect->PE=(int *)ckalloc(GAPARRSIZE*sizeof(int));//newly added
		//fprintf(stderr,"creating array for PEs in a connection.\n");
		int i;
		for(i=0;i<weight;i++){
			//connect->PE[i]=gap;
			//fprintf(stderr,"inputting a PE with estimated gap size %d\n",gap);
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

static int in1PE(unsigned int e1,int pre_pos,unsigned int bal_e2,int pos,int insert_size)
{
	int gap,realpeSize;
	unsigned int bal_e1,e2;
	if(e1==bal_e2){
		//ignorePE1++;
		return -1;  //orientation wrong
	}

	bal_e1 = getTwinCtg(e1);
	e2 = getTwinCtg(bal_e2);
	if(e1==e2){
		realpeSize = contig_array[e1].length + overlaplen - pre_pos - pos;
		if(realpeSize>0){
			//peSUM += realpeSize;
			onsameCtgPE++;
			if((int)contig_array[e1].length>insert_size){
				int *item = (int *)stackPush(isStack);
				(*item) = realpeSize;
			}
		}
		return 2;
	}

	gap = insert_size - overlaplen + pre_pos + pos - contig_array[e1].length - contig_array[e2].length;
	//fprintf(stderr,"[%s]\t%d\t%d\tgap\t%d\t%d\t%d\t%d\n",__FUNCTION__,e1,e2,gap,contig_array[e1].bal_edge,contig_array[e2].bal_edge,insert_size);
	//if(gap<-(insert_size/10)){
	//	//ignorePE2++;
	//	return 0;
	//}
	bun1AccuConnect(e1,e2,gap,1);
	bun1AccuConnect(bal_e2,bal_e1,gap,1);

	return 1;
}

static int inputPE(FILE *fp,int peGrad,char *line)
{
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
		
	//onsameCtgPE = peSUM = 0;
	PE = pes[peGrad].insertS;
	if(strlen(line)){
			sscanf(line,"%lld %d %d",&pre_readno,&pre_contigno,&pre_pos);
			//printf("first record %d %d %d\n",pre_readno,pre_contigno,pre_pos);
			if(pre_readno<=minno)
				pre_readno = -1;
	}
	else
		pre_readno = -1;
	//ignorePE1 = ignorePE2 = ignorePE3 = ignorePE4 = ignorePE5 = 0;
	//static_flag = 1;
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
				flag = in1PE(pre_contigno,pre_pos,newIndex,pos,PE);
				if(flag==1)
					count++;
			}
			pre_readno = readno;
			pre_contigno = newIndex;
			pre_pos = pos;
	}
	printf("[%s]Finish loading all PEs in grad %d .\n",__FUNCTION__,peGrad);
	printf("[%s]Calculating estimated gap size for all connections .\n",__FUNCTION__);
	/*unsigned int i;
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
					//ignorePE4+=2;
					tmp->gapLen=sum/(weight-2);
					//fprintf(stderr,"estimating contigs' gap by removing max&min PE ,with max&min %d %d\n",
						//tmp->PE[maxid],tmp->PE[minid]);
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
						sum+=((avg-(long long int)tmp->PE[ii])*(avg-(long long int)tmp->PE[ii]));
					}
					
					double SD=(sqrt((double)sum/(weight-1)))*3;//just for fast
					sum=0;
					int num=0;
					for(ii=0;ii<weight;ii++){
						if(abs(tmp->PE[ii]-avg)<=SD){
							sum+=tmp->PE[ii];
							num++;
						}else{
							//ignorePE5++;
							counter++;
						}
					}
					if(num==0){
						//fprintf(stderr,"[%s]num=0 in removing exceed 3*SD(%.1f) avg(%lld)step",__FUNCTION__,SD,avg);
						for(ii=0;ii<weight;ii++){
							fprintf(stderr,"%d\t",tmp->PE[ii]);
						}
					}
					tmp->gapLen=sum/num;
					//fprintf(stderr,"estimating contigs' gap by removing PE exceeding 3*SD ,removing %d PEs\n",counter);
				}else if(tmp->weightNotInherit<=2){
					int weight=tmp->weightNotInherit;
					int sum=0;
					int ii;
					for(ii=0;ii<weight;ii++){
						sum+=tmp->PE[ii];
					}
					tmp->gapLen=sum/weight;
					//fprintf(stderr,"weight too small , directly estimate gap size.\n");
				//}
				//fprintf(stderr,"finish %d connection.\n",i);
				//free((void *)tmp->PE);
				tmp=tmp->next;
			}
	}*/
	//printf("%d PEs with insert size %d attached, %d + %d + %d ignored\n",count,PE,ignorePE1,ignorePE2,ignorePE3);
	fprintf(stderr,"[%s]%d PEs of insert size %d loaded .\n",__FUNCTION__,count,PE);
	//fprintf(stderr,"[%s]PEs discarded:%d because of wrong orientation,%d too close,%d too far,\n",__FUNCTION__,ignorePE1,ignorePE2,ignorePE3);
	//fprintf(stderr,"[%s]%d deleted by removing max&min , %d not fall in 3*SD.\n",__FUNCTION__,ignorePE4,ignorePE5);
	//printf("[%s]%d PEs of insert size %d loaded .\n",__FUNCTION__,count,PE);
	//printf("[%s]PEs discarded :%d because of wrong orientation,%d too close,%d too far ,\n",__FUNCTION__,ignorePE1,ignorePE2,ignorePE3);
	//printf("[%s]%d deleted by removing max&min , %d not fall in 3*SD .\n",__FUNCTION__,ignorePE4,ignorePE5);
	
	/*if(onsameCtgPE>0){
		//printf("estimated PE size %lli, by %d pairs\n",peSUM/onsameCtgPE,onsameCtgPE);
		int SD=0;
		int avg=calcuIS(isStack,&SD);
		printf("[%s]%d PE attached on same contig with estimated insert size %d SD %d .\n",__FUNCTION__,onsameCtgPE,avg,SD);
	}*/
	//printf("on contigs longer than %d, %d pairs found,",PE,isStack->item_c);
	//printf("insert_size estimated: %d\n",calcuIS(isStack));
	//freeStack(isStack);
	return count;
}

int call_bundle(){
	char name[256],*line;
	FILE *fp,*linkF;
	int i;
	int flag=0;
	unsigned int j;

	loadUpdatedEdges(graphfile);
	
	//sprintf(name,"%s.bundle",graphfile);

	linkF = ckopen(name,"w");

	if(!pes)
		loadPEgrads(graphfile);

	sprintf(name,"%s.readOnContig",graphfile);
	fp = ckopen(name,"r");
	
	lineLen = 1024;
	line = (char *)ckalloc(lineLen*sizeof(char));

	fgets(line,lineLen,fp);
	line[0] = '\0';

	//printf("\n");
	newCntCounter = 0;
	//createCntMemManager();
	//createCntLookupTable();
	/*int *length_array = (unsigned int *)ckalloc((num_ctg+1)*sizeof(unsigned int));
	//use length_array to change info in index_array
	for(i=1;i<=num_ctg;i++)
		length_array[i] = 0;

	for(i=1;i<=num_ctg;i++){
		if(index_array[i]>0)
			length_array[index_array[i]] = i;
	}
	for(i=1;i<=num_ctg;i++)
		index_array[i] = length_array[i];
	*/
	for(i=0;i<gradsCounter;i++){

		createCntMemManager();
		createCntLookupTable();
		//
		flag += inputPE(fp,i,line);
		//sprintf(name,"%d.bundle",i);

		//printf("%lld new connections\n",newCntCounter/2);
		/*if(!flag){
			destroyConnectMem();
			deleteCntLookupTable();
			for(j=1;j<=num_ctg;j++)
				contig_array[j].downwardConnect = NULL;
			//printf("\n");
			continue;
		}*/
		flag = 0;
		//linkF= ckopen(name,"w");
		//outputBundle(linkF, pes[i].insertS);
		
		for(j=1;j<=num_ctg;j++){
			CONNECT *tmp=contig_array[j].downwardConnect;
			while(tmp){
				free((void *)tmp->PE);
				tmp=tmp->next;
			}
			contig_array[j].downwardConnect = NULL;
		}
		//destroyConnectMem();
		//deleteCntLookupTable();

		fclose(linkF);
	}
	 
	outputBundle(linkF,1);
	destroyConnectMem();
	deleteCntLookupTable();

	free((void *)line);
	fclose(fp);
	//fclose(linkF);
	printf("[%s]all PEs attached\n",__FUNCTION__);

	return 0;
}

void outputBundle(FILE *fp, int insertS)
{
	unsigned int i,bal_ctg,bal_toCtg;
	CONNECT *cnts,*temp_cnt;
	//printf("outputLinks, %d contigs\n",num_ctg);
	for(i=1;i<=num_ctg;i++){
		cnts = contig_array[i].downwardConnect;
		bal_ctg = getTwinCtg(i);
		//fprintf(stderr,"contig %d.\n",i);
		while(cnts){
			if(cnts->weightNotInherit<=bund_threshold){
				cnts = cnts->next;
				continue;
			}
			//fprintf(stderr,"with contig %d.\n",cnts->contigID);
			//fprintf(fp,"%-10d %-10d\t%d\t%d\t%d\n"
			//,i,cnts->contigID,cnts->gapLen,cnts->weight,insertS);
			/*int st1,st2,ed1,ed2,len1,len2,gap;
			len1=contig_array[i].length+overlaplen;
			len2=contig_array[cnts->contigID].length+overlaplen;
			gap=-cnts->gapLen;
			if(len1<gap){
				st1=0;
				if(len2<gap){
					ed1=len1-gap+len2;
					st2=gap-len1;
					ed2=len2;
				}else{
					ed1=len1;
					st2=gap-len1;
					ed2=gap;
			}
			}else{
				st1=len1+overlaplen-gap;
				st2=0;
				if(len2<gap){
					ed1=st1+len2;
					ed2=len2;
				}else{
					ed2=gap;
					ed1=len1;
				}
			}
			
			unsigned int id1,id2;
			id1=index_array[i];
			id2=index_array[cnts->contigID];*/
			/*if((id1/2+1)==1194){
				int ii;
				fprintf(stdout,"\n");
				for(ii=0;ii<cnts->weightNotInherit;++ii){
					fprintf(stdout,"%d ",cnts->PE[ii]);
				}
				fprintf(stdout,"\n");
			}*/
			/*if(isSmallerThanTwin(id1)){
				if(isSmallerThanTwin(id2)){
					fprintf(fp,"%u\t%d\t%u\t%d\t%d\n",id1/2+1,len1,id2/2+1,len2,cnts->gapLen,cnts->weightNotInherit);
				
				}else{
					fprintf(fp,"%u\t%d\t%u\t%d\t%d\n",id1/2+1,len1,-id2/2,len2,cnts->gapLen,cnts->weightNotInherit);
				}
			}else{
				if(isSmallerThanTwin(id2)){
					fprintf(fp,"%u\t%d\t%u\t%d\t%d\n",-id1/2,len1,id2/2+1,len2,cnts->gapLen,cnts->weightNotInherit);
				}else{
					fprintf(fp,"%u\t%d\t%u\t%d\t%d\n",-id1/2,len1,-id2/2,len2,cnts->gapLen,cnts->weightNotInherit);
				}
			}*/
			//int ii=0;
			//int weight=cnts->weightNotInherit;
			//for(;ii<weight;++ii){
			//	fprintf(fp,"%d\t%d\t%d\t",icnts->gapLen);
			//}
			if(cnts->gapLen<0){
				fprintf(fp,"%d\t%d\t%d\n",i,cnts->contigID,cnts->gapLen);
			}

			//fprintf(fp,"\n");
			cnts->weightNotInherit = 0;
			
			bal_toCtg = getTwinCtg(cnts->contigID);
			temp_cnt = getCntBetween(bal_toCtg,bal_ctg);
			if(temp_cnt)
				temp_cnt->weightNotInherit = 0;
			
			cnts = cnts->next;
		}
	}
}


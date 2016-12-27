#include "stdinc.h" 
#include "newhash.h"
#include "extfunc.h" 
#include "extvab.h" 
#include "dfibHeap.h"
#include "fibHeap.h"
#include "darray.h"

#define CNBLOCKSIZE 10000
#define MAXC 10000
#define MAXCinBetween 200

#define MaxNodeInSub 10000
#define GapLowerBound -2000
#define GapUpperBound 300000

//static boolean static_f=0;


static int gapCounter;
static int orienCounter;
static int throughCounter;

static DARRAY *solidArray;
static DARRAY *tempArray;

static int solidCounter;

static CTGinHEAP ctg4heapArray[MaxNodeInSub+1];  // index in this array are put to heaps, start from 1
static unsigned int nodesInSub[MaxNodeInSub];
static int nodeDistance[MaxNodeInSub];
static int nodeCounter;

static unsigned int nodesInSubInOrder[MaxNodeInSub];
static int nodeDistanceInOrder[MaxNodeInSub];

static DARRAY *scaf3,*scaf5;
static DARRAY *gap3,*gap5;

static unsigned int downstreamCTG[MAXCinBetween];
static unsigned int upstreamCTG[MAXCinBetween];
static int dsCtgCounter;
static int usCtgCounter;

static CONNECT *checkConnect(unsigned int from_c,unsigned int to_c);
static int maskPuzzle(int num_connect,unsigned int contigLen);
static void freezing();
static boolean checkOverlapInBetween(double tolerance);
static int setConnectDelete(unsigned int from_c,unsigned int to_c,char flag,boolean cleanBinding);
static int setConnectWP(unsigned int from_c,unsigned int to_c,char flag);

static void general_linearization(boolean strict);
static void debugging2();
static void smallScaf();
static void detectBreakScaf();
static boolean checkSimple(DARRAY *ctgArray,int count);
static void checkCircle();

//find the only connection involved in connection binding
static CONNECT *getBindCnt(unsigned int ctg)
{
	CONNECT *ite_cnt;
	CONNECT *bindCnt=NULL;
	CONNECT *temp_cnt=NULL;
	CONNECT *temp3_cnt=NULL;
	int count = 0;
	int count2 = 0;
	int count3 = 0;
	
	ite_cnt = contig_array[ctg].downwardConnect;
	while(ite_cnt){
		if(ite_cnt->nextInScaf){
			count++;
			bindCnt = ite_cnt;
		}
		if(ite_cnt->prevInScaf){
			temp_cnt = ite_cnt;
			count2++;
		}
		if(ite_cnt->singleInScaf){
			temp3_cnt = ite_cnt;
			count3++;
		}	
		ite_cnt = ite_cnt->next;	
	}
	if(count==1)
		return bindCnt;
	
	if(count==0&&count2==1)
		return temp_cnt;
	if(count==0&&count2==0&&count3==1)
		return temp3_cnt;
	return NULL;
}

static void createAnalogousCnt(unsigned int sourceStart,
		CONNECT *originCnt, int gap,
			unsigned int targetStart,unsigned int targetStop)
{
	CONNECT *temp_cnt;
	unsigned int balTargetStart=getTwinCtg(targetStart);
	unsigned int balTargetStop=getTwinCtg(targetStop);

	unsigned int balSourceStart = getTwinCtg(sourceStart);
	unsigned int balSourceStop = getTwinCtg(originCnt->contigID);
	
	originCnt->deleted = 1;
	temp_cnt = getCntBetween(balSourceStop,balSourceStart);
	temp_cnt->deleted = 1;

	if(gap<GapLowerBound){
		gapCounter++;
		return;
	}
	temp_cnt = add1Connect(targetStart,targetStop,gap,originCnt->weight,1);
	if(temp_cnt)
		temp_cnt->inherit = 1;
	temp_cnt = add1Connect(balTargetStop,balTargetStart,gap,originCnt->weight,1);
	if(temp_cnt)
		temp_cnt->inherit = 1;
}
// increase #long_pe_support for a conncet by 1
static void add1LongPEcov(unsigned int fromCtg,unsigned int toCtg,int weight)
{
	//check if they are on the same scaff
	if(contig_array[fromCtg].from_vt!=contig_array[toCtg].from_vt ||
		contig_array[fromCtg].to_vt!=contig_array[toCtg].to_vt){
		printf("Warning from add1LongPEcov: contig %d and %d not on the same scaffold\n",
			fromCtg,toCtg);
		return;
	}
	if(contig_array[fromCtg].indexInScaf>=contig_array[toCtg].indexInScaf){
		printf("Warning from add1LongPEcov: wrong about order between contig %d and %d\n",
			fromCtg,toCtg);
		return;
	}
	CONNECT *bindCnt;
	unsigned int prevCtg = fromCtg;
	bindCnt = getBindCnt(fromCtg);
	while(bindCnt){
		if(bindCnt->maxGap + weight<=1000)
			bindCnt->maxGap += weight;
		else
			bindCnt->maxGap = 1000;

		if(fromCtg==0&&toCtg==0)
			printf("link (%d %d ) covered by link (%d %d), wt %d\n",
				prevCtg,bindCnt->contigID,fromCtg,toCtg,weight);
		if(bindCnt->contigID==toCtg)
			break;
		prevCtg = bindCnt->contigID;
		bindCnt = bindCnt->nextInScaf;
	}
	unsigned int bal_fc = getTwinCtg(fromCtg);
	unsigned int bal_tc = getTwinCtg(toCtg);
	bindCnt = getBindCnt(bal_tc);
	prevCtg = bal_tc;
	while(bindCnt){
		if(bindCnt->maxGap + weight<=1000)
			bindCnt->maxGap += weight;
		else
			bindCnt->maxGap = 1000;
		if(fromCtg==0&&toCtg==0)
			printf("link (%d %d ) covered by link (%d %d), wt %d\n",
				prevCtg,bindCnt->contigID,fromCtg,toCtg,weight);
		if(bindCnt->contigID==bal_fc)
			return;
		prevCtg = bindCnt->contigID;
		bindCnt = bindCnt->nextInScaf;
	}
	printf("Warning from add1LongPEcov: not reach the end (%d %d) (B)\n",bal_tc,bal_fc);
}

// for long pair ends, move the connections along scaffolds established by shorter pair ends till reach the ends
static void downSlide()
{
	fprintf(stderr,"[%s]entering this function.\n",__FUNCTION__);
	int len=0,gap;
	unsigned int i;
	CONNECT *ite_cnt,*bindCnt,*temp_cnt;
	unsigned int bottomCtg,topCtg,bal_i;
	unsigned int targetCtg,bal_target;
	boolean getThrough,orienConflict;
	int slideLen,slideLen2;

	orienCounter = throughCounter = 0;
	for(i=1;i<=num_ctg;i++){
		if(contig_array[i].mask||!contig_array[i].downwardConnect)
			continue;
		bindCnt = getBindCnt(i);
		if(!bindCnt)
			continue;
		bal_i = getTwinCtg(i);
		len = slideLen = 0;
		bottomCtg = i;	

		//find the last unmasked contig in this binding
		while(bindCnt->nextInScaf){
			len += bindCnt->gapLen + contig_array[bindCnt->contigID].length;
			if(contig_array[bindCnt->contigID].mask==0){
				bottomCtg = bindCnt->contigID;
				slideLen = len;
			}
			bindCnt = bindCnt->nextInScaf;
		}
		len += bindCnt->gapLen + contig_array[bindCnt->contigID].length;

		if(contig_array[bindCnt->contigID].mask==0||bottomCtg==0){
			bottomCtg = bindCnt->contigID;
			slideLen = len;
		}
		//check each connetion from long pair ends
		ite_cnt = contig_array[i].downwardConnect;
		while(ite_cnt){
			if(ite_cnt->deleted||ite_cnt->mask||ite_cnt->singleInScaf
				||ite_cnt->nextInScaf||ite_cnt->prevInScaf||ite_cnt->inherit){
				ite_cnt = ite_cnt->next;	
				continue;
			}
			targetCtg = ite_cnt->contigID;
			if(contig_array[i].from_vt==contig_array[targetCtg].from_vt){  // on the same scaff
				if(contig_array[i].indexInScaf>contig_array[targetCtg].indexInScaf)
					orienCounter++;
				else
					throughCounter++;
			
				setConnectDelete(i,ite_cnt->contigID,1,0);
				ite_cnt = ite_cnt->next;	
				continue;

			}
			//check if this connection conflicts with previous scaffold orientationally 
			temp_cnt = getBindCnt(targetCtg);
			orienConflict = 0;
			if(temp_cnt){
				while(temp_cnt->nextInScaf){
					if(temp_cnt->contigID==i){
						orienConflict = 1;
						printf("Warning from downSlide: still on the same scaff: %d and %d\n"
							,i,targetCtg);
						printf("on scaff %d and %d\n",
							contig_array[i].from_vt,contig_array[targetCtg].from_vt);
						printf("on bal_scaff %d and %d\n",
							contig_array[bal_target].to_vt,contig_array[bal_i].to_vt);
						break;
					}
					temp_cnt = temp_cnt->nextInScaf;
				}
				if(temp_cnt->contigID==i)
					orienConflict = 1;
			}
			if(orienConflict){
				orienCounter++;
				setConnectDelete(i,ite_cnt->contigID,1,0);
				ite_cnt = ite_cnt->next;	
				continue;
			}
			//find the most top contig along previous scaffold starting with the target contig of this connection
			bal_target = getTwinCtg(targetCtg);
			slideLen2 = 0;
			if(contig_array[targetCtg].mask==0){
				topCtg = bal_target;
			}else{
				topCtg = 0;
			}
		
			temp_cnt = getBindCnt(bal_target);
			getThrough = len = 0;
			if(temp_cnt){
				//find the last contig in this binding
				while(temp_cnt->nextInScaf){
					//check if this route reaches bal_i
					if(temp_cnt->contigID==bal_i){
						printf("Warning from downSlide: (B) still on the same scaff: %d and %d (%d and %d)\n",
								i,targetCtg,bal_target,bal_i);
						printf("on scaff %d and %d\n",
							contig_array[i].from_vt,contig_array[targetCtg].from_vt);
						printf("on bal_scaff %d and %d\n",
							contig_array[bal_target].to_vt,contig_array[bal_i].to_vt);
						getThrough = 1;
						break;
					}
					len += temp_cnt->gapLen + contig_array[temp_cnt->contigID].length;
					if(contig_array[temp_cnt->contigID].mask==0){
						topCtg = temp_cnt->contigID;
						slideLen2 = len;
					}
					temp_cnt = temp_cnt->nextInScaf;
				}	
				len += temp_cnt->gapLen + contig_array[temp_cnt->contigID].length;
				if(contig_array[temp_cnt->contigID].mask==0||topCtg==0){
					topCtg = temp_cnt->contigID;
					slideLen2 = len;
				}
				if(temp_cnt->contigID==bal_i)
					getThrough = 1;
				else 
					topCtg = getTwinCtg(topCtg);
			}else
				topCtg = targetCtg;

			if(getThrough){
				throughCounter++;
				setConnectDelete(i,ite_cnt->contigID,1,0);
				ite_cnt = ite_cnt->next;	
				continue;
			}
			//add a connection between bottomCtg and topCtg
			gap = ite_cnt->gapLen - slideLen - slideLen2;
			if(bottomCtg!=topCtg&&!(i==bottomCtg&&targetCtg==topCtg)){
				createAnalogousCnt(i,ite_cnt,gap,bottomCtg,topCtg);
				if(contig_array[bottomCtg].mask||contig_array[topCtg].mask)
					printf("downSlide to masked contig\n");
			}
			ite_cnt = ite_cnt->next;	
		} //for each connect
	}  // for each contig
	//printf("downSliding is done...orienConflict %d, fall inside %d\n",
		//		orienCounter,throughCounter);
}

static boolean setNextInScaf(CONNECT *cnt, CONNECT *nextCnt)
{
	if(!cnt){
		printf("setNextInScaf: empty pointer\n");
		return 0;
	}
	if(!nextCnt){
		cnt->nextInScaf = nextCnt;	
		return 1;
	}
	if(cnt->mask||cnt->deleted){
		printf("setNextInScaf: cnt is masked or deleted\n");
		return 0;
	}
	if(nextCnt->deleted||nextCnt->mask){
		printf("setNextInScaf: nextCnt is masked or deleted\n");
		return 0;
	}
	cnt->nextInScaf = nextCnt;
	return 1;
}

static boolean setPrevInScaf(CONNECT *cnt, boolean flag)
{
	if(!cnt){
		printf("setPrevInScaf: empty pointer\n");
		return 0;
	}
	if(!flag){
		cnt->prevInScaf = flag;
		return 1;
	}
	if(cnt->mask||cnt->deleted){
		printf("setPrevInScaf: cnt is masked or deleted\n");
		return 0;
	}
	cnt->prevInScaf = flag;
	return 1;
}

/*
connect A is upstream to B, replace A with C
from_c
 			> branch_c - to_c
from_c_new
*/
static void substitueUSinScaf(CONNECT *origin, unsigned int from_c_new)
{
	if(!origin||!origin->nextInScaf)
		return;

	unsigned int branch_c, to_c;
	unsigned int bal_branch_c, bal_to_c;
	unsigned int bal_from_c_new = getTwinCtg(from_c_new);
	CONNECT *bal_origin,*bal_nextCNT,*prevCNT,*bal_prevCNT;
	
	
	branch_c = origin->contigID;
	to_c = origin->nextInScaf->contigID;
	bal_branch_c = getTwinCtg(branch_c);
	bal_to_c = getTwinCtg(to_c);
	
	prevCNT = checkConnect(from_c_new,branch_c);
	bal_nextCNT = checkConnect(bal_to_c,bal_branch_c);
	if(!bal_nextCNT){
		printf("substitueUSinScaf: no connect between %d and %d\n",bal_to_c,bal_branch_c);	
		return;
	}
	bal_origin = bal_nextCNT->nextInScaf;
	bal_prevCNT = checkConnect(bal_branch_c,bal_from_c_new);

	setPrevInScaf(bal_nextCNT->nextInScaf,0);
	setNextInScaf(prevCNT,origin->nextInScaf);
	setNextInScaf(bal_nextCNT,bal_prevCNT);
	setPrevInScaf(bal_prevCNT,1);

	setNextInScaf(origin,NULL);
	setPrevInScaf(bal_origin,0);
}

/*
connect B is downstream to C, replace B with A
                       to_c
from_c - branch_c < 
				       to_c_new
*/
static void substitueDSinScaf(CONNECT *origin, unsigned int branch_c, unsigned int to_c_new)
{
	if(!origin||!origin->prevInScaf)
		return;

	unsigned int to_c;
	unsigned int bal_branch_c, bal_to_c,bal_to_c_new;
	unsigned int from_c,bal_from_c;
	CONNECT *bal_origin,*prevCNT,*bal_prevCNT;
	CONNECT *nextCNT,*bal_nextCNT;
	
	
	to_c = origin->contigID;
	bal_branch_c = getTwinCtg(branch_c);
	bal_to_c = getTwinCtg(to_c);
	bal_origin = getCntBetween(bal_to_c,bal_branch_c);
	if(!bal_origin){
		printf("substitueDSinScaf: no connect between %d and %d\n",bal_to_c,bal_branch_c);	
		return;
	}
	bal_from_c = bal_origin->nextInScaf->contigID;
	from_c = getTwinCtg(bal_from_c);
	bal_to_c_new = getTwinCtg(to_c_new);
	
	prevCNT = checkConnect(from_c,branch_c);
	nextCNT = checkConnect(branch_c,to_c_new);
	setNextInScaf(prevCNT,nextCNT);
	setPrevInScaf(nextCNT,1);

	bal_nextCNT = checkConnect(bal_to_c_new,bal_branch_c);
	bal_prevCNT = checkConnect(bal_branch_c,bal_from_c);

	setNextInScaf(bal_nextCNT,bal_prevCNT);
	setPrevInScaf(origin,0);
	setNextInScaf(bal_origin,NULL);
}

static int validConnect(unsigned int ctg, CONNECT *preCNT)
{
	if(preCNT&&preCNT->nextInScaf)
		return 1;

	CONNECT *cn_temp;
	int count=0;
	if(!contig_array[ctg].downwardConnect)	
		return count;
	cn_temp = contig_array[ctg].downwardConnect;
	while(cn_temp){
		if(!cn_temp->deleted&&!cn_temp->mask)
			count++;
		cn_temp = cn_temp->next;
	}
	return count;
}

static CONNECT *getNextContig(unsigned int ctg, CONNECT *preCNT, boolean *exception)
{
	CONNECT *cn_temp,*retCNT=NULL;
	int count=0,valid_in;
	unsigned int nextCtg,bal_ctg;
	
	*exception = 0;
	if(preCNT&&preCNT->nextInScaf){
		if(preCNT->contigID!=ctg)
			printf("pre cnt does not lead to %d\n",ctg);
		nextCtg = preCNT->nextInScaf->contigID;
		cn_temp = getCntBetween(ctg,nextCtg);
		if(cn_temp&&(cn_temp->mask||cn_temp->deleted)){
			printf("getNextContig: arc(%d %d) twin (%d %d) with mask %d deleted %d\n"
							,ctg,nextCtg,getTwinCtg(nextCtg),getTwinCtg(ctg)
							,cn_temp->mask,cn_temp->deleted);
			if(!cn_temp->prevInScaf)
				printf("not even has a prevInScaf\n");
			cn_temp = getCntBetween(getTwinCtg(nextCtg),
							getTwinCtg(ctg));
			if(!cn_temp->nextInScaf)
				printf("its twin cnt not  has a nextInScaf\n");
			fflush(stdout);
			*exception = 1;
		}else
			return preCNT->nextInScaf;
	}

	bal_ctg = getTwinCtg(ctg);
	valid_in = validConnect(bal_ctg,NULL);
	if(valid_in>1)
		return NULL;
	if(!contig_array[ctg].downwardConnect)
		return NULL;
	cn_temp = contig_array[ctg].downwardConnect;
	while(cn_temp){
		if(cn_temp->mask||cn_temp->deleted){
			cn_temp = cn_temp->next;
			continue;
		}
		count++;
		if(count==1)
			retCNT = cn_temp;
		else if(count==2)
			return NULL;
		cn_temp = cn_temp->next;
	}
	return retCNT;
}

// get the valid connect between 2 given ctgs
static CONNECT *checkConnect(unsigned int from_c,unsigned int to_c)
{
	CONNECT *cn_temp=getCntBetween(from_c,to_c);
	if(!cn_temp)
		return NULL;
	if(!cn_temp->mask&&!cn_temp->deleted)
		return cn_temp;
	return NULL;
}

static int setConnectMask(unsigned int from_c,unsigned int to_c,char mask)
{
	CONNECT *cn_temp,*cn_bal,*cn_ds,*cn_us;
	unsigned int bal_fc = getTwinCtg(from_c);
	unsigned int bal_tc = getTwinCtg(to_c);
	unsigned int ctg3,bal_ctg3;

	cn_temp = getCntBetween(from_c,to_c);
	cn_bal = getCntBetween(bal_tc,bal_fc);
	if(!cn_temp||!cn_bal){
		return 0;
	}
	cn_temp->mask = mask;
	cn_bal->mask = mask;
	if(!mask)
		return 1;

	if(cn_temp->nextInScaf){ //undo the binding
		setPrevInScaf(cn_temp->nextInScaf,0);
		ctg3 = cn_temp->nextInScaf->contigID;
		setNextInScaf(cn_temp,NULL);
		bal_ctg3 = getTwinCtg(ctg3);
		cn_ds = getCntBetween(bal_ctg3,bal_tc);
		setNextInScaf(cn_ds,NULL);
		setPrevInScaf(cn_bal,0);
	}

	// ctg3 -> from_c -> to_c
	// bal_ctg3 <- bal_fc <- bal_tc
	if(cn_bal->nextInScaf){
		setPrevInScaf(cn_bal->nextInScaf,0);
		bal_ctg3 = cn_bal->nextInScaf->contigID;
		setNextInScaf(cn_bal,NULL);
		ctg3 = getTwinCtg(bal_ctg3);
		cn_us = getCntBetween(ctg3,from_c);
		setNextInScaf(cn_us,NULL);
		setPrevInScaf(cn_temp,0);
	}

	return 1;
}


static boolean setConnectUsed(unsigned int from_c,unsigned int to_c,char flag)
{
	CONNECT *cn_temp,*cn_bal;
	unsigned int bal_fc = getTwinCtg(from_c);
	unsigned int bal_tc = getTwinCtg(to_c);

	cn_temp = getCntBetween(from_c,to_c);
	cn_bal = getCntBetween(bal_tc,bal_fc);
	if(!cn_temp||!cn_bal){
		return 0;
	}
	cn_temp->used = flag;
	cn_bal->used = flag;
	
	return 1;
}

static int setConnectWP(unsigned int from_c,unsigned int to_c,char flag)
{
	CONNECT *cn_temp,*cn_bal;
	unsigned int bal_fc = getTwinCtg(from_c);
	unsigned int bal_tc = getTwinCtg(to_c);

	cn_temp = getCntBetween(from_c,to_c);
	cn_bal = getCntBetween(bal_tc,bal_fc);
	if(!cn_temp||!cn_bal){
		return 0;
	}
	cn_temp->weakPoint = flag;
	cn_bal->weakPoint = flag;
	//fprintf(stderr,"contig %d and %d, weakPoint %d\n",from_c,to_c,cn_temp->weakPoint);
	//fprintf(stderr,"contig %d and %d, weakPoint %d\n",bal_tc,bal_fc,cn_bal->weakPoint);
	return 1;
}

static int setConnectDelete(unsigned int from_c,unsigned int to_c,char flag,boolean cleanBinding)
{
	CONNECT *cn_temp,*cn_bal;
	unsigned int bal_fc = getTwinCtg(from_c);
	unsigned int bal_tc = getTwinCtg(to_c);

	cn_temp = getCntBetween(from_c,to_c);
	cn_bal = getCntBetween(bal_tc,bal_fc);

	if(!cn_temp||!cn_bal){
		return 0;
	}
	cn_temp->deleted = flag;
	cn_bal->deleted = flag;
	if(!flag)
		return 1;
	if(cleanBinding){
		cn_temp->prevInScaf = 0;
		cn_temp->nextInScaf = NULL;
		cn_bal->prevInScaf = 0;
		cn_bal->nextInScaf = NULL;
	}
	return 1;
}

static void maskContig(unsigned int ctg,boolean flag)
{
	unsigned int bal_ctg,ctg2,bal_ctg2;
	CONNECT *cn_temp;

	bal_ctg = getTwinCtg(ctg);
	cn_temp = contig_array[ctg].downwardConnect;
	while(cn_temp){
		if(cn_temp->mask||cn_temp->prevInScaf||cn_temp->nextInScaf||cn_temp->singleInScaf){
			cn_temp = cn_temp->next;
			continue;
		}
		ctg2 = cn_temp->contigID;
		setConnectMask(ctg,ctg2,flag);
		cn_temp = cn_temp->next;
	}
	// bal_ctg2 <- bal_ctg
	cn_temp = contig_array[bal_ctg].downwardConnect;
	while(cn_temp){
		if(cn_temp->mask||cn_temp->prevInScaf||cn_temp->nextInScaf||cn_temp->singleInScaf){
			cn_temp = cn_temp->next;
			continue;
		}
		bal_ctg2 = cn_temp->contigID;
		setConnectMask(bal_ctg,bal_ctg2,flag);
		cn_temp = cn_temp->next;
	}
	
	contig_array[ctg].mask = flag;
	contig_array[bal_ctg].mask = flag;
}

static int maskPuzzle(int num_connect,unsigned int contigLen)
{
	int in_num,out_num,flag=0,puzzleCounter=0;
	unsigned int i,bal_i;

	for(i=1;i<=num_ctg;i++){
		if(contigLen&&contig_array[i].length>contigLen)
			break;
		if(contig_array[i].mask)
			continue;
		bal_i = getTwinCtg(i);
		in_num = validConnect(bal_i,NULL);
		out_num = validConnect(i,NULL);
		if((in_num>1||out_num>1)&&(in_num+out_num>=num_connect)){
			flag++;
			maskContig(i,1);	
		}
		in_num = validConnect(bal_i,NULL);
		out_num = validConnect(i,NULL);
		if(in_num>1||out_num>1){
			puzzleCounter++;
			//debugging2(i);
		}
	
		if(isSmallerThanTwin(i))
			i++;
	}
	//printf("Masked %d contigs, %d puzzle left\n",flag,puzzleCounter);
	return flag;
}

static void deleteWeakCnt(int cut_off)
{
	unsigned int i;
	CONNECT *cn_temp1;
	int weaks=0,counter=0;
	//fprintf(stderr,"[%s]entering this function. num_ctg=%d\n",__FUNCTION__,num_ctg);
	for(i=1;i<=num_ctg;i++){
		//fprintf(stderr,"[%s]iterating %d.\n",__FUNCTION__,i);
		cn_temp1 = contig_array[i].downwardConnect;
		while(cn_temp1){
			if(!cn_temp1->mask&&!cn_temp1->deleted&&!cn_temp1->nextInScaf
				&&!cn_temp1->singleInScaf&&!cn_temp1->prevInScaf){
				counter++;
			}
			if(cn_temp1->weak&&cn_temp1->deleted&&cn_temp1->weight>=cut_off){
				cn_temp1->deleted = 0;
				cn_temp1->weak = 0;
			}
			else if(!cn_temp1->deleted&&cn_temp1->weight>0&&cn_temp1->weight<cut_off
					&&!cn_temp1->nextInScaf&&!cn_temp1->prevInScaf){
				cn_temp1->deleted = 1;
				cn_temp1->weak = 1;
				if(cn_temp1->singleInScaf)
					cn_temp1->singleInScaf = 0;
				if(!cn_temp1->mask)
					weaks++;
			}
			cn_temp1 = cn_temp1->next;
		}

	}
	fprintf(stderr,"[%s]%d connects doesn't meet weight threshold .\n",__FUNCTION__,weaks);
	checkCircle();
}

//check if one contig is linearly connected to the other ->C1->C2...
static int linearC2C(unsigned int starter,CONNECT *cnt2c1,unsigned int c2,int min_dis,int max_dis)
{
	int out_num,in_num;
	CONNECT *prevCNT,*cnt,*cn_temp;
	unsigned int c1,bal_c1,ctg,bal_c2;
	int len=0;
	unsigned int bal_start = getTwinCtg(starter);
	boolean excep;

	c1 = cnt2c1->contigID;

	if(c1==c2){
		printf("linearC2C: c1(%d) and c2(%d) are the same contig\n",c1,c2);
		return -1;
	}

	bal_c1 = getTwinCtg(c1);
	in_num = validConnect(bal_c1,NULL);
	if(in_num>1)
		return 0;

	dsCtgCounter = 1;
	usCtgCounter = 0;
	downstreamCTG[dsCtgCounter++] = c1;
	bal_c2 = getTwinCtg(c2);
	upstreamCTG[usCtgCounter++] = bal_c2;
	// check if c1 is linearly connected to c2 by pe connections
	cnt = prevCNT = cnt2c1;
	while((cnt=getNextContig(c1,prevCNT,&excep))!=NULL){
		c1 = cnt->contigID;
		len += cnt->gapLen+contig_array[c1].length;
		if(c1==c2)
			return 1;
		
		if(len>max_dis||c1==starter||c1==bal_start)
			return 0;
		downstreamCTG[dsCtgCounter++] = c1;
		if(dsCtgCounter>=MAXCinBetween){
			printf("%d downstream contigs, start at %d, max_dis %d, current dis %d\n"
				,dsCtgCounter,starter,max_dis,len);
			return 0;
		}
		prevCNT = cnt;
	}
	out_num = validConnect(c1,NULL);
	if(out_num)
		return 0;
	

	//find the most upstream contig to c2
	cnt = prevCNT = NULL; 
	ctg = bal_c2;
	while((cnt=getNextContig(ctg,prevCNT,&excep))!=NULL){
		ctg = cnt->contigID;
		len += cnt->gapLen+contig_array[ctg].length;
		if(len>max_dis||ctg==starter||ctg==bal_start)
			return 0;

		prevCNT = cnt;
		upstreamCTG[usCtgCounter++] = ctg;
		if(usCtgCounter>=MAXCinBetween){
			printf("%d upstream contigs, start at %d, max_dis %d, current dis %d\n"
				,usCtgCounter,starter,max_dis,len);
			return 0;
		}
	}
	if(dsCtgCounter+usCtgCounter>MAXCinBetween){
		printf("%d downstream and %d upstream contigs\n",dsCtgCounter,usCtgCounter);
		return 0;
	}
	out_num = validConnect(ctg,NULL);
	if(out_num){
		return 0;
	}
		
	c2 = getTwinCtg(ctg);
	min_dis -= len;
	max_dis -= len;  
	if(c1==c2||c1==ctg||max_dis<0)
		return 0;
		
	cn_temp = getCntBetween(c1,c2);
	if(cn_temp){
		setConnectMask(c1,c2,0);
		setConnectDelete(c1,c2,0,0);
		return 1;
	}
	len = (min_dis+max_dis)/2 >= 0 ? (min_dis+max_dis)/2 : 0;
	cn_temp = allocateCN(c2,len);
	if(cntLookupTable)
		putCnt2LookupTable(c1,cn_temp);
	cn_temp->weight = 0;  // special connect from the original graph
	cn_temp->next = contig_array[c1].downwardConnect;
	contig_array[c1].downwardConnect = cn_temp;
		
	bal_c1 = getTwinCtg(c1);
	bal_c2 = getTwinCtg(c2);

	cn_temp = allocateCN(bal_c1,len);
	if(cntLookupTable)
		putCnt2LookupTable(bal_c2,cn_temp);
	cn_temp->weight = 0;  // special connect from the original graph
	cn_temp->next = contig_array[bal_c2].downwardConnect;
	contig_array[bal_c2].downwardConnect = cn_temp;
	return 1;
}
//catenate upstream contig array and downstream contig array to solidArray
static void catUsDsContig()
{
	int i;
	
	for(i=0;i<dsCtgCounter;i++)
		*(unsigned int *)darrayPut(solidArray,i) = downstreamCTG[i];

	for(i=usCtgCounter-1;i>=0;i--){
		*(unsigned int *)darrayPut(solidArray,dsCtgCounter++) = getTwinCtg(upstreamCTG[i]); 
	}

	solidCounter = dsCtgCounter;
}

//binding the connections between contigs in solidArray 
static void consolidate()
{
	int i,j;
	CONNECT *prevCNT=NULL;
	CONNECT *cnt;
	unsigned int to_ctg;
	unsigned int from_ctg = *(unsigned int *)darrayGet(solidArray,0);

	for(i=1;i<solidCounter;i++){
		to_ctg = *(unsigned int *)darrayGet(solidArray,i);
		cnt = checkConnect(from_ctg,to_ctg);
		if(!cnt){
			printf("consolidate A: no connect from %d to %d\n",
					from_ctg,to_ctg);
			for(j=0;j<solidCounter;j++)
				printf("%d-->",*(unsigned int *)darrayGet(solidArray,j));
			printf("\n");
			return;
		}
		cnt->singleInScaf = solidCounter==2 ? 1:0;
		if(prevCNT){
			setNextInScaf(prevCNT,cnt);
			setPrevInScaf(cnt,1);
		}
		prevCNT = cnt;
		from_ctg = to_ctg;
	}

	//the reverse complementary path
	from_ctg = getTwinCtg(*(unsigned int*)darrayGet(solidArray,solidCounter-1));
	prevCNT = NULL;
	for(i=solidCounter-2;i>=0;i--){
		to_ctg = getTwinCtg(*(unsigned int *)darrayGet(solidArray,i)); 
		cnt = checkConnect(from_ctg,to_ctg);
		if(!cnt){
			printf("consolidate B: no connect from %d to %d\n",from_ctg,to_ctg);
			return;
		}
		cnt->singleInScaf = solidCounter==2 ? 1:0;
		if(prevCNT){
			setNextInScaf(prevCNT,cnt);
			setPrevInScaf(cnt,1);
		}
		prevCNT = cnt;
		from_ctg = to_ctg;
	}
	
}

static void debugging1(unsigned int ctg1,unsigned int ctg2)
{
	CONNECT *cn1;
	cn1 = getCntBetween(ctg1,ctg2);
	if(cn1){
		printf("(%d,%d) mask %d deleted %d w %d,singleInScaf %d\n",
				ctg1,ctg2,cn1->mask,cn1->deleted,cn1->weight,cn1->singleInScaf);
		if(cn1->nextInScaf)
			printf("%d->%d->%d\n",ctg1,ctg2,cn1->nextInScaf->contigID);
		if(cn1->prevInScaf)
			printf("*->%d->%d\n",ctg1,ctg2);
		else if(!cn1->nextInScaf)	
			printf("NULL->%d->%d->NULL\n",ctg1,ctg2);
	}else
		printf("%d -X- %d\n",ctg1,ctg2);
}
//remove transitive connections which cross linear paths (these paths may be broken)
//if a->b->c and a->c, mask a->c
static void removeTransitive()
{
	unsigned int i,bal_ctg;
	int flag=1,out_num,in_num,count,min,max,linear;
	CONNECT *cn_temp,*cn1=NULL,*cn2=NULL;

	while(flag){
		flag = 0;
		for(i=1;i<=num_ctg;i++){
			if(contig_array[i].mask)
				continue;
			out_num = validConnect(i,NULL);
			if(out_num!=2)
				continue;
			cn_temp = contig_array[i].downwardConnect;
			count = 0;
			while(cn_temp){
				if(cn_temp->deleted||cn_temp->mask){
					cn_temp = cn_temp->next;
					continue;
				}
				count++;
				if(count==1)
					cn1 = cn_temp;
				else if(count==2){
					cn2 = cn_temp;
				}else // count > 2
					break;
				
				cn_temp = cn_temp->next;
			}	
			if(count>2){
				printf("%d valid connections from ctg %d\n",count,i);
				continue;
			}
			if(cn1->gapLen>cn2->gapLen){
				cn_temp = cn1;
				cn1 = cn2;
				cn2 = cn_temp;
			}  //make sure cn1 is closer to contig i than cn2
			if(cn1->prevInScaf&&cn2->prevInScaf)
				continue;
			bal_ctg = getTwinCtg(cn2->contigID);
			in_num = validConnect(bal_ctg,NULL);
			if(in_num>2)
				continue;
			min = cn2->gapLen - cn1->gapLen - contig_array[cn1->contigID].length - ins_size_var/2;
			max = cn2->gapLen - cn1->gapLen - contig_array[cn1->contigID].length + ins_size_var/2;
			
			if(max<0)
				continue;
			//temprarily delete cn2
			setConnectDelete(i,cn2->contigID,1,0);
			linear = linearC2C(i,cn1,cn2->contigID,min,max);
			if(linear!=1){
				setConnectDelete(i,cn2->contigID,0,0);
				continue;
			}else{
				downstreamCTG[0] = i;
				catUsDsContig();
				if(!checkSimple(solidArray,solidCounter))
					continue;
				cn1 = getCntBetween(*(unsigned int *)darrayGet(solidArray,solidCounter-2),
						*(unsigned int *)darrayGet(solidArray,solidCounter-1));
				if(cn1&&cn1->nextInScaf&&cn2->nextInScaf){
					setConnectDelete(i,cn2->contigID,0,0);
					continue;
				}
				consolidate();
				if(cn2->prevInScaf)
					substitueDSinScaf(cn2,*(unsigned int *)darrayGet(solidArray,0),
						*(unsigned int *)darrayGet(solidArray,1));
				if(cn2->nextInScaf)
					substitueUSinScaf(cn2,*(unsigned int *)darrayGet(solidArray,solidCounter-2));
				flag++;
			}
		}  //for each contig
		//printf("a remove transitive lag, %d connections removed\n",flag);
	}

}

//get repeat contigs back into the scaffold according to connected unique contigs on both sides
/*
	A  ------  D 
       > [i] <
	B          E 
*/
static void debugging2(unsigned int ctg)
{
	CONNECT *cn1 = contig_array[ctg].downwardConnect;
	while(cn1){
		if(cn1->nextInScaf)
			fprintf(stderr,"with nextInScaf,");
		if(cn1->prevInScaf)
			fprintf(stderr,"with prevInScaf,");
		fprintf(stderr,"%u >> %d, mask %d deleted %d, inherit %d, singleInScaf %d\n",
			ctg,cn1->contigID,cn1->mask,cn1->deleted,cn1->inherit,cn1->singleInScaf);
		cn1 = cn1->next;	
	}
}
static void debugging()
{
/*
	debugging1(1777,1468);
	debugging2(8065);
	debugging2(8066);
*/
}

static void simplifyCnt()
{
	removeTransitive();
	debugging();
	general_linearization(1);
	debugging();
}

static int getIndexInArray(unsigned int node)
{
	int index;
	for(index=0;index<nodeCounter;index++)
		if(nodesInSub[index]==node)
			return index;
	return -1;
}

static boolean putNodeIntoSubgraph(FibHeap *heap, int distance, unsigned int node, int index)
{

	int pos = getIndexInArray(node);
	if(pos>0){
		//printf("exists\n");
		return 0;
	}
	if(index>=MaxNodeInSub)
		return -1;
	insertNodeIntoHeap(heap,distance,node);
	nodesInSub[index] = node;
	nodeDistance[index] = distance;
	return 1;
}

static boolean putChainIntoSubgraph(FibHeap *heap,int distance,unsigned int node,int *index,CONNECT *prevC)
{
	unsigned int ctg = node;
	CONNECT *nextCnt;
	boolean excep,flag;
	int counter = *index;

	while(1){
		nextCnt=getNextContig(ctg,prevC,&excep);
		if(excep||!nextCnt){
			*index = counter;
			return 1;
		}
		ctg = nextCnt->contigID;
		distance += nextCnt->gapLen + ctg;
		flag = putNodeIntoSubgraph(heap,distance,ctg,counter);
		if(flag<0)
			return 0;
		if(flag>0)
			counter++;
		prevC = nextCnt;
	}
}
//check if nodes in subgraph have a potential heter form
static boolean check_het_overlap(double tolerance){

	int i,gap,overlap_point;
	unsigned int node;
	int len_sum,over3_len,over3_sum;
	boolean flag=0;
	len_sum=0;
	over3_len=0;
	over3_sum=0;
	for(i=1;i<=nodeCounter;i++){
		node = ctg4heapArray[i].ctgID;
		len_sum += contig_array[node].length;
	}
	if(len_sum<1)
		return 2;
	for(i=1;i<nodeCounter;i++){
		gap = ctg4heapArray[i+1].dis - ctg4heapArray[i].dis
				- contig_array[ctg4heapArray[i+1].ctgID].length;
		if(gap>0){
			flag=0;
		}
		else{
			if(flag){
				over3_len=ctg4heapArray[i+1].dis - overlap_point
					- contig_array[ctg4heapArray[i+1].ctgID].length;
				over3_sum+=over3_len;
				if((double)over3_sum/len_sum>tolerance)
					return 0;
			}
			flag=1;
			overlap_point=ctg4heapArray[i].dis;
		}
	}

	return 2;
}

// check if a contig is unique by trying to line its downstream/upstream nodes together
static boolean checkUnique(unsigned int node,double tolerance)
{
	CONNECT *ite_cnt;
	unsigned int currNode;
	int distance;
	int popCounter = 0;
	boolean flag;

	currNode = node;
	FibHeap *heap = newFibHeap();

	putNodeIntoSubgraph(heap, 0, currNode, 0);
	nodeCounter = 1;
	ite_cnt = contig_array[currNode].downwardConnect;
	while(ite_cnt){
		if(ite_cnt->deleted||ite_cnt->mask){
			ite_cnt = ite_cnt->next;	
			continue;
		}
		currNode = ite_cnt->contigID;
		distance = ite_cnt->gapLen + contig_array[currNode].length;
		flag = putNodeIntoSubgraph(heap, distance, currNode, nodeCounter);
		if(flag<0){
			destroyHeap(heap);	
			return 0;
		}
		if(flag>0)
			nodeCounter++;

		flag = putChainIntoSubgraph(heap,distance,currNode,&nodeCounter,ite_cnt);
		if(!flag){
			destroyHeap(heap);	
			return 0;
		}

		ite_cnt = ite_cnt->next;	
	}
	if(nodeCounter<=2){  // no more than 2 valid connections
		destroyHeap(heap);	
		return 1;
	}
	
	while((currNode=removeNextNodeFromHeap(heap))!=0)
		nodesInSubInOrder[popCounter++] = currNode;
		
	destroyHeap(heap);

	flag = checkOverlapInBetween(tolerance);
	if(flag==1){
		return 1;
	}else{
		flag = check_het_overlap(0.02);//check the heter form 
	}
	return flag;
}

//find longest path and break the other
static void process_ds_contig(unsigned int ctg){
	unsigned int target=ctg4heapArray[nodeCounter].ctgID;
	//int boarder = ctg4heapArray[nodeCounter].dis;
	boolean excep;
	CONNECT *route=contig_array[ctg].downwardConnect;
	CONNECT *max_route=route;
	
	int max_dis=0;
	
	boolean end_flag=0;
	while(route){
		
		int dis=0;
		CONNECT *tmp_cnt=route;	
		while(tmp_cnt){
			dis+=route->gapLen+contig_array[route->contigID].length;
			if(route->contigID==target){
				end_flag=1;
				break;
			}
			tmp_cnt=getNextContig(route->contigID,tmp_cnt,&excep);
		}
		if(dis>max_dis){
			max_dis=dis;
			max_route=route;
		}
		if(end_flag){
			max_route=route;
			break;
		}
		route=route->next;
	}
	//delete connect except max_route
	route=contig_array[ctg].downwardConnect;
	while(route){
		if(route!=max_route){
			setConnectMask(ctg,route->contigID,1);
		}
		route=route->next;
	}
	
}
static void process_us_contig(unsigned int ctg){
	unsigned int target=ctg4heapArray[1].ctgID;
	//int boarder = ctg4heapArray[1].dis;
	boolean excep;
	CONNECT *route=contig_array[ctg].downwardConnect;
	CONNECT *min_route=route;
	
	int min_dis=0;
	
	boolean end_flag=0;
	while(route){
		
		int dis=0;
		CONNECT *tmp_cnt=route;	
		while(tmp_cnt){
			dis-=route->gapLen+contig_array[route->contigID].length;
			if(route->contigID==target){
				end_flag=1;
				break;
			}
			tmp_cnt=getNextContig(route->contigID,tmp_cnt,&excep);
		}
		if(dis<min_dis){
			min_dis=dis;
			min_route=route;
		}
		if(end_flag){
			min_route=route;
			break;
		}
		route=route->next;
	}
	//delete connect except min_route
	route=contig_array[ctg].downwardConnect;
	while(route){
		if(route!=min_route){
			setConnectMask(ctg,route->contigID,1);
		}
		route=route->next;
	}
	
}

//mask contigs with downstream and/or upstream can not be lined 
static void maskRepeat()
{
	int in_num,out_num,flagA,flagB;
	int counter = 0;
	int puzzleCounter = 0;
	unsigned int i,bal_i;
	int het_counter = 0;
	for(i=1;i<=num_ctg;i++){
		if(contig_array[i].mask)
			continue;
		bal_i = getTwinCtg(i);
		in_num = validConnect(bal_i,NULL);
		out_num = validConnect(i,NULL);
		if(in_num>1||out_num>1)
			puzzleCounter++;
		else{
			if(isSmallerThanTwin(i))
				i++;
			continue;
			
		}

		if(contig_array[i].cvg>2*cvgAvg){
			counter++;
			maskContig(i,1);
			//printf("thick mask contig %d and %d\n",i,bal_i);
			if(isSmallerThanTwin(i))
				i++;
			continue;
		}

		if(in_num>1)
			flagA = checkUnique(bal_i,OverlapPercent);
		else
			flagA = 1;
		if(out_num>1)
			flagB = checkUnique(i,OverlapPercent);
		else
			flagB = 1;

		if(flagA==0||flagB==0){
			counter++;
			maskContig(i,1);
		}else{
			if(flagA==2){//us find longest path
				process_us_contig(bal_i);
			}
			if(flagB==2){//ds find longest path
				process_ds_contig(i);
			}	
		}
		if(flagA==2||flagB==2)
			het_counter++;

		if(isSmallerThanTwin(i))
			i++;
	}
	printf("[%s]%d contigs masked from %d puzzles\n",__FUNCTION__,counter,puzzleCounter);
	printf("[%s]%d processed as heterozygous .\n",__FUNCTION__,het_counter);
}


static void ordering(boolean deWeak,boolean downS, boolean nonlinear, char *infile)
{
			//debugging();
			if(downS){
				downSlide();
				//debugging();
				if(deWeak)
					deleteWeakCnt(weakPE);
			}else{
				if(deWeak)
					deleteWeakCnt(weakPE);
			}
			//output_scaf(infile);
			//debugging();
			//printf("variance for insert size %d\n",ins_size_var);
			simplifyCnt();
			//debugging();

			maskRepeat();
			//debugging();
			simplifyCnt();	
			
			if(nonlinear){
				//printf("non-strict linearization\n");
				general_linearization(0);
				//linearization(0,0);
			}
			//maskRepeat();//???
			
			maskPuzzle(2,0);
			//debugging();
			freezing();
			//debugging();
	
}

//check if contigs next to each other have reasonable overlap
boolean checkOverlapInBetween(double tolerance)
{
	int i,gap;
	int index;
	unsigned int node;
	int lenSum,lenOlp;
	lenSum = lenOlp = 0;
	for(i=0;i<nodeCounter;i++){
		node = nodesInSubInOrder[i];
		lenSum += contig_array[node].length;
		index = getIndexInArray(node);
		nodeDistanceInOrder[i] = nodeDistance[index];
	}
	if(lenSum<1)
		return 1;
	for(i=0;i<nodeCounter-1;i++){
		gap = nodeDistanceInOrder[i+1] - nodeDistanceInOrder[i]
				- contig_array[nodesInSubInOrder[i+1]].length;
		if(-gap>0)
			lenOlp += -gap;
		//if(-gap>ins_size_var)
		if((double)lenOlp/lenSum>tolerance)
			return 0;
	}
	return 1;
}


/*********   the following codes are for freezing current scaffolds   ****************/
//set connections between contigs in a array to used or not
//meanwhile set mask to the opposite value
static boolean setUsed(unsigned int start,unsigned int *array,int max_steps,boolean flag)
{
	unsigned int prevCtg = start;
	unsigned int twinA,twinB;
	int j;
	CONNECT *cnt;
	boolean usedFlag=0;
	// save 'used' to 'checking'
	prevCtg = start;
	for(j=0;j<max_steps;j++){
		if(array[j]==0)
			break;
		cnt = getCntBetween(prevCtg,array[j]);
		if(!cnt){
			printf("setUsed: no connect between %d and %d\n",prevCtg,array[j]);	
			prevCtg = array[j];
			continue;
		}
		if(cnt->used==flag||cnt->nextInScaf||cnt->prevInScaf||cnt->singleInScaf){
			return 1;
		}
		cnt->checking = cnt->used;
		twinA = getTwinCtg(prevCtg);
		twinB = getTwinCtg(array[j]);
		cnt = getCntBetween(twinB,twinA);
		if(cnt)
			cnt->checking = cnt->used;
		prevCtg = array[j];
	}
	// set used to flag
	prevCtg = start;
	for(j=0;j<max_steps;j++){
		if(array[j]==0)
			break;
		cnt = getCntBetween(prevCtg,array[j]);
		if(!cnt){
			prevCtg = array[j];
			continue;
		}
		if(cnt->used==flag){
			usedFlag = 1;
			break;  
		}
		cnt->used = flag;
		twinA = getTwinCtg(prevCtg);
		twinB = getTwinCtg(array[j]);
		cnt = getCntBetween(twinB,twinA);
		if(cnt)
			cnt->used = flag;
		prevCtg = array[j];
	}
	// set mask to 'NOT flag' or set used to original value
	prevCtg = start;
	for(j=0;j<max_steps;j++){
		if(array[j]==0)
			break;
		cnt = getCntBetween(prevCtg,array[j]);
		if(!cnt){
			prevCtg = array[j];
			continue;
		}
		if(!usedFlag)
			cnt->mask = 1-flag;
		else
			cnt->used = cnt->checking;
		twinA = getTwinCtg(prevCtg);
		twinB = getTwinCtg(array[j]);
		cnt = getCntBetween(twinB,twinA);
		cnt->used = 1-flag;
		if(!usedFlag)
			cnt->mask = 1-flag;
		else
			cnt->used = cnt->checking;
		prevCtg = array[j];
	}	
	return usedFlag;
}
// break down scaffolds poorly supported by longer PE
static void recoverMask()
{
	unsigned int i,ctg,bal_ctg,start,finish;
	int num3,num5,j,t;
	CONNECT *bindCnt,*cnt;
	int min,max,max_steps=5,num_route,length;
	int tempCounter,recoverCounter=0;
	boolean multiUSE,change;

	for(i=1;i<=num_ctg;i++)
		contig_array[i].flag = 0;

	so_far = (unsigned int *)ckalloc(max_n_routes*sizeof(unsigned int));
	found_routes = (unsigned int **)ckalloc(max_n_routes*sizeof(unsigned int *));
	for(j=0;j<max_n_routes;j++)
		found_routes[j] = (unsigned int *)ckalloc(max_steps*sizeof(unsigned int));
	for(i=1;i<=num_ctg;i++){
		if(contig_array[i].flag||contig_array[i].mask||!contig_array[i].downwardConnect)
			continue;
		bindCnt = getBindCnt(i);
		if(!bindCnt)
			continue;
		//first scan get the average coverage by longer pe
		num5 = num3 = 0;
		ctg = i;
		*(unsigned int *)darrayPut(scaf5,num5++) = i;
		contig_array[i].flag = 1;
		contig_array[getTwinCtg(i)].flag = 1;
		while(bindCnt){
			if(bindCnt->used)
				break;
			setConnectUsed(ctg,bindCnt->contigID,1);
			ctg = bindCnt->contigID;
			*(unsigned int *)darrayPut(scaf5,num5++) = ctg;
			bal_ctg = getTwinCtg(ctg);
			contig_array[ctg].flag = 1;
			contig_array[bal_ctg].flag = 1;
			bindCnt = bindCnt->nextInScaf;
		}

		ctg = getTwinCtg(i);
		bindCnt = getBindCnt(ctg);
		while(bindCnt){
			if(bindCnt->used)
				break;
			setConnectUsed(ctg,bindCnt->contigID,1);
			ctg = bindCnt->contigID;
			bal_ctg = getTwinCtg(ctg);
			contig_array[ctg].flag = 1;
			contig_array[bal_ctg].flag = 1;
			*(unsigned int *)darrayPut(scaf3,num3++) = bal_ctg;
			bindCnt = bindCnt->nextInScaf;
		}
		if(num5+num3<2)
			continue;
		tempCounter = solidCounter = 0;
		for(j=num3-1;j>=0;j--)
			*(unsigned int *)darrayPut(tempArray,tempCounter++) = 	
						*(unsigned int *)darrayGet(scaf3,j);
		for(j=0;j<num5;j++)
			*(unsigned int *)darrayPut(tempArray,tempCounter++) = 
					*(unsigned int *)darrayGet(scaf5,j);
		change = 0;
		for(t=0;t<tempCounter-1;t++){
			*(unsigned int *)darrayPut(solidArray,solidCounter++) = 
					*(unsigned int *)darrayGet(tempArray,t);
			start = *(unsigned int *)darrayGet(tempArray,t);
			finish = *(unsigned int *)darrayGet(tempArray,t+1);
			num_route = num_trace = 0;
			cnt = checkConnect(start,finish);
			if(!cnt){
				printf("Warning from recoverMask: no connection (%d %d), start at %d\n",
					start,finish,i);	
				cnt = getCntBetween(start,finish);
				if(cnt)
					debugging1(start,finish);
				continue;
			}
			length = cnt->gapLen + contig_array[finish].length;
			min = length - 1.5*ins_size_var;
			max = length + 1.5*ins_size_var;
			traceAlongMaskedCnt(finish,start,max_steps,min,max,0,0,&num_route);
			if(finish==start){
				for(j=0;j<tempCounter;j++)
					printf("->%d",*(unsigned int *)darrayGet(tempArray,j));
				printf(": start at %d\n",i);
			}
			
			if(num_route==1){
				for(j=0;j<max_steps;j++)
					if(found_routes[0][j]==0)
						break;
				if(j<1)
					continue;
				//check if connects have been used more than once
				multiUSE = setUsed(start,found_routes[0],max_steps,1);
				if(multiUSE)
					continue;
				for(j=0;j<max_steps;j++){
					if(j+1==max_steps||found_routes[0][j+1]==0)
						break;
					*(unsigned int *)darrayPut(solidArray,solidCounter++) = found_routes[0][j];
					contig_array[found_routes[0][j]].flag = 1;
					contig_array[getTwinCtg(found_routes[0][j])].flag = 1;
				}
				recoverCounter += j;
				setConnectDelete(start,finish,1,1);
				change = 1;
			}  //end if num_route=1
		}  // for each gap
		*(unsigned int *)darrayPut(solidArray,solidCounter++) = 
				*(unsigned int *)darrayGet(tempArray,tempCounter-1);
		if(change)
			consolidate();
	}

	//printf("%d contigs recovered\n",recoverCounter);
	fflush(stdout);

	for(i=1;i<=num_ctg;i++){
		cnt = contig_array[i].downwardConnect;
		while(cnt){
			cnt->used = 0;
			cnt->checking = 0;
			cnt = cnt->next;	
		}
	}	

	for(j=0;j<max_n_routes;j++)
		free((void *)found_routes[j]);
	free((void *)found_routes);
	free((void *)so_far);
}


// A -> B -> C -> D  un-bind link B->C to link A->B and B->C
// A' <- B' <- C' <- D'
static void unBindLink(unsigned int CB,unsigned int CC)
{
	//fprintf(stderr,"Unbind link (%d %d) to others...\n",CB,CC);
	CONNECT *cnt1 = getCntBetween(CB,CC);
	if(!cnt1)
		return;
	if(cnt1->singleInScaf)
		cnt1->singleInScaf = 0;
	CONNECT *cnt2 = getCntBetween(getTwinCtg(CC),getTwinCtg(CB));
	if(!cnt2)
		return;
	if(cnt2->singleInScaf)
		cnt2->singleInScaf = 0;
	if(cnt1->nextInScaf){
		unsigned int CD = cnt1->nextInScaf->contigID;
		cnt1->nextInScaf->prevInScaf = 0;
		cnt1->nextInScaf = NULL;
		CONNECT *cnt3 = getCntBetween(getTwinCtg(CD),getTwinCtg(CC));
		if(cnt3)
			cnt3->nextInScaf = NULL;
		cnt2->prevInScaf = 0;
	}
	if(cnt2->nextInScaf){
		unsigned int bal_CA = cnt2->nextInScaf->contigID;
		cnt2->nextInScaf->prevInScaf = 0;
		cnt2->nextInScaf = NULL;
		CONNECT *cnt4 = getCntBetween(getTwinCtg(bal_CA),CB);
		if(cnt4)
			cnt4->nextInScaf = NULL;
		cnt1->prevInScaf = 0;
	}
}

static void freezing()
{
	int num5,num3;
	unsigned int ctg,bal_ctg;
	unsigned int i;
	int j,t;
	CONNECT *cnt,*prevCNT,*nextCnt;
	boolean excep;

	for(i=1;i<=num_ctg;i++){
		contig_array[i].flag = 0;
		contig_array[i].from_vt = 0;
		contig_array[i].to_vt = 0;
		cnt = contig_array[i].downwardConnect;
		while(cnt){
			cnt->used = 0;
			cnt->checking = 0;
			cnt->singleInScaf = 0;
			cnt = cnt->next;	
		}
	}
	
	for(i=1;i<=num_ctg;i++){
		if(contig_array[i].flag||contig_array[i].mask)
			continue;

		if(!contig_array[i].downwardConnect||!validConnect(i,NULL)){
			continue;
		}
		
		num5 = num3 = 0;
		ctg = i;
		*(unsigned int *)darrayPut(scaf5,num5++) = i;
		contig_array[i].flag = 1;
		contig_array[getTwinCtg(i)].flag = 1;
		prevCNT = NULL;
		cnt = getNextContig(ctg,prevCNT,&excep);
		while(cnt){
			if(contig_array[cnt->contigID].flag){
				unBindLink(ctg,cnt->contigID);	
				break;
			}
			nextCnt=getNextContig(cnt->contigID,cnt,&excep);
			setConnectUsed(ctg,cnt->contigID,1);
			ctg = cnt->contigID;
			*(unsigned int *)darrayPut(scaf5,num5++) = ctg;
			bal_ctg = getTwinCtg(ctg);
			contig_array[ctg].flag = 1;
			contig_array[bal_ctg].flag = 1;
			prevCNT = cnt;
			cnt = nextCnt;
		}

		ctg = getTwinCtg(i);
		if(num5>=2)
			prevCNT = checkConnect(getTwinCtg(*(unsigned int *)darrayGet(scaf5,1)),ctg);
		else
			prevCNT = NULL;
		cnt = getNextContig(ctg,prevCNT,&excep);
		while(cnt){
			if(contig_array[cnt->contigID].flag){
				unBindLink(ctg,cnt->contigID);	
				break;
			}
			nextCnt=getNextContig(cnt->contigID,cnt,&excep);
			setConnectUsed(ctg,cnt->contigID,1);
			ctg = cnt->contigID;
			bal_ctg = getTwinCtg(ctg);
			contig_array[ctg].flag = 1;
			contig_array[bal_ctg].flag = 1;
			*(unsigned int *)darrayPut(scaf3,num3++) = bal_ctg;
			prevCNT = cnt;
			cnt = nextCnt;
		}
		if(num5+num3<2)
			continue;
		solidCounter = 0;
		for(j=num3-1;j>=0;j--)
			*(unsigned int *)darrayPut(solidArray,solidCounter++) = 	
						*(unsigned int *)darrayGet(scaf3,j);
		for(j=0;j<num5;j++)
			*(unsigned int *)darrayPut(solidArray,solidCounter++) = 
					*(unsigned int *)darrayGet(scaf5,j);
		unsigned int firstCtg = 0;
		unsigned int lastCtg = 0;
		unsigned int firstTwin = 0;
		unsigned int lastTwin = 0;
		for(t=0;t<solidCounter;t++)
			if(!contig_array[*(unsigned int *)darrayGet(solidArray,t)].mask){
				firstCtg = *(unsigned int *)darrayGet(solidArray,t);
				break;
			}
		for(t=solidCounter-1;t>=0;t--)
			if(!contig_array[*(unsigned int *)darrayGet(solidArray,t)].mask){
				lastCtg = *(unsigned int *)darrayGet(solidArray,t);
				break;
			}
		if(firstCtg==0||lastCtg==0){
			printf("scaffold start at %d, stop at %d, freezing began with %d\n",firstCtg,lastCtg,i);
			for(j=0;j<solidCounter;j++)
				printf("->%d(%d %d)",*(unsigned int *)darrayGet(solidArray,j)
					,contig_array[*(unsigned int *)darrayGet(solidArray,j)].mask
					,contig_array[*(unsigned int *)darrayGet(solidArray,j)].flag);
			printf("\n");
		}else{
			firstTwin = getTwinCtg(firstCtg);
			lastTwin = getTwinCtg(lastCtg);
		}
		for(t=0;t<solidCounter;t++){
			unsigned int ctg = *(unsigned int *)darrayGet(solidArray,t);
			if(contig_array[ctg].from_vt>0){
				contig_array[ctg].mask = 1;
				contig_array[getTwinCtg(ctg)].mask = 1;
				printf("Repeat: contig %d (%d) appears more than once\n",ctg,getTwinCtg(ctg));
			}else{
				contig_array[ctg].from_vt = firstCtg;
				contig_array[ctg].to_vt = lastCtg;
				contig_array[ctg].indexInScaf = t+1;
				contig_array[getTwinCtg(ctg)].from_vt = lastTwin;
				contig_array[getTwinCtg(ctg)].to_vt = firstTwin;
				contig_array[getTwinCtg(ctg)].indexInScaf = solidCounter-t;
			}
		}
		consolidate();
	}

	//printf("Freezing is done....\n");
	fflush(stdout);

	for(i=1;i<=num_ctg;i++){
		if(contig_array[i].flag)
			contig_array[i].flag = 0;
		
		if(contig_array[i].from_vt==0){
			contig_array[i].from_vt = i;
			contig_array[i].to_vt = i;
		}
		cnt = contig_array[i].downwardConnect;
		while(cnt){
			cnt->used = 0;
			cnt->checking = 0;
			cnt = cnt->next;	
		}
	}	

}

/************** codes below this line are for pulling the scaffolds out ************/
void output1gap(FILE *fo,int max_steps)
{
	int i,len,seg;
	len = seg = 0;

	for(i=0;i<max_steps-1;i++){
		if(found_routes[0][i+1]==0)
			break;
		len += contig_array[found_routes[0][i]].length;
		seg++;
	}
	fprintf(fo,"GAP %d %d",len,seg);
	for(i=0;i<max_steps-1;i++){
		if(found_routes[0][i+1]==0)
			break;
		
		fprintf(fo," %d",found_routes[0][i]);
	}
	fprintf(fo,"\n");
}

static int weakCounter;

static boolean printCnts(FILE *fp,unsigned int ctg)
{
	CONNECT *cnt = contig_array[ctg].downwardConnect;
	boolean flag = 0,ret=0;
	unsigned int bal_ctg = getTwinCtg(ctg);
	unsigned int linkCtg;
	if(isSameAsTwin(ctg))
		return ret;
	CONNECT *bindCnt = getBindCnt(ctg);
	if(bindCnt&&bindCnt->bySmall&&bindCnt->weakPoint){
		weakCounter++;
		fprintf(fp,"\tWP");
		ret = 1;
	}

	while(cnt){
		if(cnt->weight&&!cnt->inherit){
			if(!flag){
				flag = 1;
				fprintf(fp,"\t#DOWN ");
			}
			linkCtg = cnt->contigID;
			if(isLargerThanTwin(linkCtg))
				linkCtg = getTwinCtg(linkCtg);
	
			fprintf(fp,"%d:%d:%d ",index_array[linkCtg],cnt->weight,cnt->gapLen);
		}
		cnt = cnt->next;
	}
	flag = 0;
	cnt = contig_array[bal_ctg].downwardConnect;
	while(cnt){
		if(cnt->weight&&!cnt->inherit){
			if(!flag){
				flag = 1;
				fprintf(fp,"\t#UP ");
			}
			linkCtg = cnt->contigID;
			if(isLargerThanTwin(linkCtg))
				linkCtg = getTwinCtg(linkCtg);
	
			fprintf(fp,"%d:%d:%d ",index_array[linkCtg],cnt->weight,cnt->gapLen);
		}
		cnt = cnt->next;
	}
	fprintf(fp,"\n");
	return ret;
}

void scaffolding(unsigned int len_cut,char *outfile)
{
	unsigned int prev_ctg,ctg,bal_ctg,*length_array,count=0,num_lctg=0;
	unsigned int i,max_steps=5;
	int num5,num3,j,len,flag,num_route,gap_c=0;
	short gap=0;
	long long sum=0,N50,N90;
	FILE *fp,*fo=NULL;
	char name[256];
	CONNECT *cnt,*prevCNT,*nextCnt;
	boolean excep,weak;
	weakCounter = 0;

	so_far = (unsigned int *)ckalloc(max_n_routes*sizeof(unsigned int));
	found_routes = (unsigned int **)ckalloc(max_n_routes*sizeof(unsigned int*));
	for(j=0;j<max_n_routes;j++)
		found_routes[j] = (unsigned int*)ckalloc(max_steps*sizeof(unsigned int));
	
	length_array = (unsigned int *)ckalloc((num_ctg+1)*sizeof(unsigned int));
	//use length_array to change info in index_array
	for(i=1;i<=num_ctg;i++)
		length_array[i] = 0;

	for(i=1;i<=num_ctg;i++){
		if(index_array[i]>0)
			length_array[index_array[i]] = i;
	}
	for(i=1;i<=num_ctg;i++)
		index_array[i] = length_array[i];  //contig i with original index: index_array[i]

	orig2new = 0;

	sprintf(name,"%s.scaf",outfile);
	fp = ckopen(name,"w");
	sprintf(name,"%s.scaf_gap",outfile);
	fo = ckopen(name,"w");

	scaf3 = (DARRAY *)createDarray(1000,sizeof(unsigned int));
	scaf5 = (DARRAY *)createDarray(1000,sizeof(unsigned int));
	gap3 = (DARRAY *)createDarray(1000,sizeof(int));
	gap5 = (DARRAY *)createDarray(1000,sizeof(int));

	for(i=1;i<=num_ctg;i++)
		contig_array[i].flag = 0;
	for(i=1;i<=num_ctg;i++){
		if(contig_array[i].length+(unsigned int)overlaplen>=len_cut)
			num_lctg++;
		else
			continue;
		if(contig_array[i].flag||contig_array[i].mask||!contig_array[i].downwardConnect||!validConnect(i,NULL))
			continue;
		
		num5 = num3 = 0;
		ctg = i;
		//printf("%d",i);
		*(unsigned int *)darrayPut(scaf5,num5++) = i;
		contig_array[i].flag = 1;
		bal_ctg = getTwinCtg(ctg);
		contig_array[bal_ctg].flag = 1;
		len = contig_array[i].length;
		prevCNT = NULL;
		cnt = getNextContig(ctg,prevCNT,&excep);
		while(cnt){
			nextCnt = getNextContig(cnt->contigID,cnt,&excep);
			if(excep&&prevCNT)
				printf("scaffolding: exception --- prev cnt from %u\n",prevCNT->contigID);
			if(nextCnt&&nextCnt->used)
				break;
			setConnectUsed(ctg,cnt->contigID,1);
			*(int *)darrayPut(gap5,num5-1) = cnt->gapLen;
			ctg = cnt->contigID;
			*(unsigned int *)darrayPut(scaf5,num5++) = ctg;
			len += cnt->gapLen+contig_array[ctg].length;
			bal_ctg = getTwinCtg(ctg);
			contig_array[ctg].flag = 1;
			contig_array[bal_ctg].flag = 1;
			prevCNT = cnt;
			cnt = nextCnt;
			//printf("->%d",ctg);
		}
		//printf("\n");

		ctg = getTwinCtg(i);
		if(num5>=2)
			prevCNT = checkConnect(getTwinCtg(*(unsigned int *)darrayGet(scaf5,1)),ctg);
		else
			prevCNT = NULL;
		//printf("%d",i);
		//fflush(stdout);
		cnt = getNextContig(ctg,prevCNT,&excep);
		while(cnt){
			nextCnt=getNextContig(cnt->contigID,cnt,&excep);
			if(excep&&prevCNT)
				printf("scaffolding: exception -- prev cnt from %u\n",prevCNT->contigID);
			if(nextCnt&&nextCnt->used)
				break;
			setConnectUsed(ctg,cnt->contigID,1);
			ctg = cnt->contigID;
			len += cnt->gapLen+contig_array[ctg].length;
			bal_ctg = getTwinCtg(ctg);
			contig_array[ctg].flag = 1;
			contig_array[bal_ctg].flag = 1;
			//printf("<-%d",bal_ctg);
			//fflush(stdout);
			*(int *)darrayPut(gap3,num3) = cnt->gapLen;
			*(unsigned int *)darrayPut(scaf3,num3++) = bal_ctg;
			prevCNT = cnt;
			cnt = nextCnt;
		}
		//printf("\n");
		len += overlaplen;
		sum += len;
		length_array[count++] = len;
		if(num5+num3<1){
			//printf("no scaffold created for contig %d\n",i);
			continue;
		}
		fprintf(fp,">scaffold%d %d %d\n",count,num5+num3,len);
		fprintf(fo,">scaffold%d %d %d\n",count,num5+num3,len);
		len = prev_ctg = 0;
		for(j=num3-1;j>=0;j--){
			if(!isLargerThanTwin(*(unsigned int *)darrayGet(scaf3,j))){
				fprintf(fp,"%-10d %-10d +   %d "
					,index_array[*(unsigned int *)darrayGet(scaf3,j)],len,
					contig_array[*(unsigned int *)darrayGet(scaf3,j)].length+overlaplen);
				weak = printCnts(fp,*(unsigned int *)darrayGet(scaf3,j));
				/*
				if(weak)
					fprintf(stderr,"scaffold%d\n",count);
				*/
			}else{
				fprintf(fp,"%-10d %-10d -   %d "
					,index_array[getTwinCtg(*(unsigned int *)darrayGet(scaf3,j))],len
					,contig_array[*(unsigned int *)darrayGet(scaf3,j)].length+overlaplen);
				weak = printCnts(fp,*(unsigned int *)darrayGet(scaf3,j));
				/*
				if(weak)
					fprintf(stderr,"scaffold%d\n",count);
				*/
			}
			if(prev_ctg){
				num_route = num_trace = 0;
				traceAlongArc(*(unsigned int *)darrayGet(scaf3,j),prev_ctg,max_steps
					,gap-ins_size_var,gap+ins_size_var,0,0,&num_route);
				if(num_route==1){
					output1gap(fo,max_steps);
					gap_c++;
				}
			}
			fprintf(fo,"%-10d %-10d\n",*(unsigned int *)darrayGet(scaf3,j),len);
			len += contig_array[*(unsigned int *)darrayGet(scaf3,j)].length + *(int *)darrayGet(gap3,j);
			prev_ctg = *(unsigned int *)darrayGet(scaf3,j);
			gap = *(int *)darrayGet(gap3,j)>0 ? *(int *)darrayGet(gap3,j):0;
		}
		for(j=0;j<num5;j++){
			if(!isLargerThanTwin(*(unsigned int *)darrayGet(scaf5,j))){
				fprintf(fp,"%-10d %-10d +   %d "
					,index_array[*(unsigned int *)darrayGet(scaf5,j)],len
					,contig_array[*(unsigned int *)darrayGet(scaf5,j)].length+overlaplen);
				weak = printCnts(fp,*(unsigned int *)darrayGet(scaf5,j));
				/*
				if(weak)
					fprintf(stderr,"scaffold%d\n",count);
				*/
			}else{
				fprintf(fp,"%-10d %-10d -   %d "
					,index_array[getTwinCtg(*(unsigned int *)darrayGet(scaf5,j))],len
					,contig_array[*(unsigned int *)darrayGet(scaf5,j)].length+overlaplen);
				weak = printCnts(fp,*(unsigned int *)darrayGet(scaf5,j));
				/*
				if(weak)
					fprintf(stderr,"scaffold%d\n",count);
				*/
			}
			if(prev_ctg){
				num_route = num_trace = 0;
				traceAlongArc(*(unsigned int *)darrayGet(scaf5,j),prev_ctg,max_steps
				,gap-ins_size_var,gap+ins_size_var,0,0,&num_route);
				if(num_route==1){
					output1gap(fo,max_steps);
					gap_c++;
				}
			}
			fprintf(fo,"%-10d %-10d\n",*(unsigned int *)darrayGet(scaf5,j),len);
			if(j<num5-1){
				len += contig_array[*(unsigned int *)darrayGet(scaf5,j)].length + 
					*(int *)darrayGet(gap5,j);
				prev_ctg = *(unsigned int *)darrayGet(scaf5,j);
				gap = *(int *)darrayGet(gap5,j)>0 ? *(int *)darrayGet(gap5,j):0;
			}
		}

	}

	freeDarray(scaf3);
	freeDarray(scaf5);
	freeDarray(gap3);
	freeDarray(gap5);

	fclose(fp);
	fclose(fo);
	//printf("\n%d scaffolds from %d contigs sum up %lldbp, with average length %lld, %d gaps filled\n"
		//	,count,num_lctg/2,sum,sum/count,gap_c);
	printf("[%s]scaffold(s) created : %d , total length : %lld.\n",__FUNCTION__,count ,sum);
	//output singleton
	for(i=1;i<=num_ctg;i++){
		if(contig_array[i].length+(unsigned int)overlaplen<len_cut||contig_array[i].flag)
			continue;
		length_array[count++] = contig_array[i].length;
		sum += contig_array[i].length;
		if(isSmallerThanTwin(i))
			i++;
	}

	// calculate N50/N90
	//printf("%d scaffolds&singleton sum up %lldbp, with average length %lld\n"
		//	,count,sum,sum/count);
	printf("[%s]total number of scaffold(s) and singleton(s) : %d, total length : %lld.\n",__FUNCTION__,count,sum);
	qsort(length_array,count,sizeof(length_array[0]),cmp_int);
	//printf("the longest is %dbp,",length_array[count-1]);
	N50 = sum*0.5;
	N90 = sum*0.9;
	sum = flag = 0;
	for(j=count-1;j>=0;j--){
		sum += length_array[j];
		if(!flag&&sum>=N50){
			printf("[%s]N50 : %d bp, ",__FUNCTION__,length_array[j]);
			flag++;
		}
		if(sum>=N90){
			printf(" N90 : %d bp\n",length_array[j]);
			break;
		}
	}
	//printf("Found %d weak points in scaffolds\n",weakCounter);
	fflush(stdout);
	free((void *)length_array);
	for(j=0;j<max_n_routes;j++)
		free((void *)found_routes[j]);
	free((void *)found_routes);
	free((void *)so_far);
}


void outputLinks(FILE *fp, int insertS)
{
	unsigned int i,bal_ctg,bal_toCtg;
	CONNECT *cnts,*temp_cnt;
	//printf("outputLinks, %d contigs\n",num_ctg);
	for(i=1;i<=num_ctg;i++){
		cnts = contig_array[i].downwardConnect;
		bal_ctg = getTwinCtg(i);
		while(cnts){
			if(cnts->weight<1){
				cnts = cnts->next;
				continue;
			}
			fprintf(fp,"%-10d %-10d\t%d\t%d\t%d\n"
				,i,cnts->contigID,cnts->gapLen,cnts->weight,insertS);
			cnts->weight = 0;

			bal_toCtg = getTwinCtg(cnts->contigID);
			temp_cnt = getCntBetween(bal_toCtg,bal_ctg);
			if(temp_cnt)
				temp_cnt->weight = 0;

			cnts = cnts->next;	
		}
	}
}

//use pe info in ascent order
void PE2Links(char *infile)
{
	fprintf(stderr,"[%s]entering this function.\n",__FUNCTION__);
	char name[256],*line;
	FILE *fp,*linkF;
	int i;
	int flag=0;
	unsigned int j;

		
	sprintf(name,"%s.links",infile);
	/*linkF = fopen(name,"r");
	if(linkF){
		printf("file %s exists, skip creating the links...\n",name);
		fclose(linkF);
		return;
	}*/

	linkF = ckopen(name,"w");

	if(!pes)
		loadPEgrads(infile);

	sprintf(name,"%s.readOnContig",infile);
	fp = ckopen(name,"r");
	
	lineLen = 1024;
	line = (char *)ckalloc(lineLen*sizeof(char));

	fgets(line,lineLen,fp);
	line[0] = '\0';

	//printf("\n");
	for(i=0;i<gradsCounter;i++){
		createCntMemManager();
		createCntLookupTable();

		newCntCounter = 0;
		flag += connectByPE_grad(fp,i,line);
		//printf("%lld new connections\n",newCntCounter/2);
		if(!flag){
			destroyConnectMem();
			deleteCntLookupTable();
			for(j=1;j<=num_ctg;j++)
				contig_array[j].downwardConnect = NULL;
			//printf("\n");
			continue;
		}
		flag = 0;
		outputLinks(linkF, pes[i].insertS);
		destroyConnectMem();
		deleteCntLookupTable();
		for(j=1;j<=num_ctg;j++)
			contig_array[j].downwardConnect = NULL;
	}
	
		
	free((void *)line);
	fclose(fp);
	fclose(linkF);
	printf("[%s]all PEs attached\n",__FUNCTION__);

}

int inputLinks(FILE *fp, int insertS,char *line)
{
	unsigned int ctg,bal_ctg,toCtg,bal_toCtg;
	int gap,wt,ins;
	unsigned int counter=0,onScafCounter=0;
	unsigned int maskCounter=0;
	if(strlen(line)){
		sscanf(line,"%d %d %d %d %d",&ctg,&toCtg,&gap,&wt,&ins);
		if(ins!=insertS)
			return counter;
		//if(contig_array[ctg].length>=ctg_short&&contig_array[toCtg].length>=ctg_short){
		if(1){
			bal_ctg = getTwinCtg(ctg);
			bal_toCtg = getTwinCtg(toCtg);
			add1Connect(ctg,toCtg,gap,wt,0);
			add1Connect(bal_toCtg,bal_ctg,gap,wt,0);
			counter++;
			if(contig_array[ctg].mask||contig_array[toCtg].mask)
				maskCounter++;
		
			if(insertS>1000&&
				contig_array[ctg].from_vt==contig_array[toCtg].from_vt&&  // on the same scaff
					contig_array[ctg].indexInScaf<contig_array[toCtg].indexInScaf){
				add1LongPEcov(ctg,toCtg,wt);
				onScafCounter++;
			}
		}
	}

	while(fgets(line,lineLen,fp)!=NULL){
		sscanf(line,"%d %d %d %d %d",&ctg,&toCtg,&gap,&wt,&ins);
		if(ins!=insertS)
		//if(ins>insertS)
			break;
		/*
		if(contig_array[ctg].length<ctg_short||contig_array[toCtg].length<ctg_short)
			continue;
		*/
		if(insertS>1000&&
			contig_array[ctg].from_vt==contig_array[toCtg].from_vt&&  // on the same scaff
				contig_array[ctg].indexInScaf<contig_array[toCtg].indexInScaf){
			add1LongPEcov(ctg,toCtg,wt);
			onScafCounter++;
		}
		bal_ctg = getTwinCtg(ctg);
		bal_toCtg = getTwinCtg(toCtg);
		add1Connect(ctg,toCtg,gap,wt,0);
		add1Connect(bal_toCtg,bal_ctg,gap,wt,0);
		counter++;
		if(contig_array[ctg].mask||contig_array[toCtg].mask)
			maskCounter++;
	}
	//printf("%d link to masked contigs, %d links on a single scaff\n",maskCounter,onScafCounter);
	return counter;
}
//use linkage info in ascent order
void Links2Scaf(char *infile)
{
	char name[256],*line;
	FILE *fp;
	int i,lib_n=0,cutoff_sum=0;
	int flag=0,flag2;
	boolean downS,nonLinear=0,smallPE=0,isPrevSmall=0,markSmall;

	if(!pes)
		loadPEgrads(infile);
	
	sprintf(name,"%s.links",infile);
	fp = ckopen(name,"r");
	
	createCntMemManager();
	createCntLookupTable();

	lineLen = 1024;
	line = (char *)ckalloc(lineLen*sizeof(char));

	fgets(line,lineLen,fp);
	line[0] = '\0';

	
	solidArray = (DARRAY *)createDarray(1000,sizeof(unsigned int));
	tempArray = (DARRAY *)createDarray(1000,sizeof(unsigned int));
	scaf3 = (DARRAY *)createDarray(1000,sizeof(unsigned int));
	scaf5 = (DARRAY *)createDarray(1000,sizeof(unsigned int));
	gap3 = (DARRAY *)createDarray(1000,sizeof(int));
	gap5 = (DARRAY *)createDarray(1000,sizeof(int));

	weakPE = 3;   //0531
	//printf("\n");
	for(i=0;i<gradsCounter;i++){
		/*if(pes[i].insertS<1000)
			isPrevSmall = 1;
		else if(pes[i].insertS>1000&&isPrevSmall){
			smallScaf();
			isPrevSmall = 0;
		}*/
		flag2 = inputLinks(fp,pes[i].insertS,line);
		//printf("Insert size %d: %d links input\n",pes[i].insertS,flag2);
		if(flag2){
			lib_n++;
			cutoff_sum += pes[i].pair_num_cut;
			weakPE=cutoff_sum;
		}
		flag += flag2;
		if(!flag){
			//printf("\n");
			continue;
		}
		if(i==gradsCounter-1|| pes[i+1].rank!=pes[i].rank){
			flag = nonLinear = downS = markSmall = 0;
	
			if(pes[i].insertS>1000&&pes[i].rank>1)
				downS = 1;
			if(pes[i].insertS<=1000)
				smallPE = 1;
			
			if(pes[i].insertS>=1000){
				ins_size_var = 50;
				//OverlapPercent = 0.05;
			}else if(pes[i].insertS>=300){
				ins_size_var = 30;
				//OverlapPercent = 0.05;
			}else{
				ins_size_var = 20;
				//OverlapPercent = 0.05;
			}
			//if(pes[i].insertS>1000)
				//weakPE = 5;
				//static_f = 1;
			//if(lib_n>0){
				//weakPE = weakPE<cutoff_sum/lib_n ? cutoff_sum/lib_n:weakPE;
				//lib_n = cutoff_sum = 0;
			//}
			
			printf("[%s]weight threshold for a connection in grad %d : %d, %d.\n",__FUNCTION__,i,weakPE,cutoff_sum);
			//printf("cut off for weight of connections : %d\n",weakPE);
			if(i==gradsCounter-1)
				nonLinear = 1;
			if(i==gradsCounter-1&&!isPrevSmall&&smallPE)
				detectBreakScaf();
			ordering(1,downS,nonLinear,infile);
			if(i==gradsCounter-1)
				recoverMask();
		}
	}
	freeDarray(tempArray);
	freeDarray(solidArray);
	freeDarray(scaf3);
	freeDarray(scaf5);
	freeDarray(gap3);
	freeDarray(gap5);

	free((void *)line);
	fclose(fp);
	//printf("all links loaded\n");

}
/* below for picking up a subgraph (with at most one node has upstream connections to the rest
									and at most one downstream connections) in general  */

// static int nodeCounter
static boolean putNodeInArray(unsigned int node, int maxNodes,int dis)
{
	if(contig_array[node].inSubGraph)
		return 1;
	int index = nodeCounter;
	if(index>maxNodes)
		return 0;
	if(contig_array[getTwinCtg(node)].inSubGraph)
		return 0;
	ctg4heapArray[index].ctgID = node;
	ctg4heapArray[index].dis = dis;
	contig_array[node].inSubGraph = 1;

	ctg4heapArray[index].ds_shut4dheap = 0;
	ctg4heapArray[index].us_shut4dheap = 0;
	ctg4heapArray[index].ds_shut4uheap = 0;
	ctg4heapArray[index].us_shut4uheap = 0;
	
	return 1;
}
	
static void setInGraph(boolean flag)
{
	int i;
	int node;
	nodeCounter = nodeCounter>MaxNodeInSub ? MaxNodeInSub:nodeCounter;
	for(i=1;i<=nodeCounter;i++){
		node = ctg4heapArray[i].ctgID;
		if(node>0)
			contig_array[node].inSubGraph = flag;
	}
}

static boolean dispatch1node(int dis,unsigned int tempNode,int maxNodes,
							FibHeap *dheap,FibHeap *uheap,int *DmaxDis,int *UmaxDis)
{
	boolean ret;
	if(dis>=0){   // put it to Dheap
		nodeCounter++;
		ret = putNodeInArray(tempNode,maxNodes,dis);
		if(!ret)
			return 0;
		insertNodeIntoHeap(dheap,dis,nodeCounter);
		if(dis>*DmaxDis)
			*DmaxDis = dis;
		return 1;
	}else{         // put it to Uheap
		nodeCounter++;
		ret = putNodeInArray(tempNode,maxNodes,dis);
		if(!ret)
			return 0;
		insertNodeIntoHeap(uheap,-dis,nodeCounter);
		int temp_len = contig_array[tempNode].length;
		if(-dis+temp_len>*UmaxDis)
			*UmaxDis = -dis+contig_array[tempNode].length;
		return -1;	
	}
	return 0;
}

static boolean canDheapWait(unsigned int currNode,int dis, int DmaxDis)
{
	if(dis<DmaxDis)
		return 0;
	else
		return 1;
}

static boolean workOnDheap(FibHeap *dheap,FibHeap *uheap,boolean *Dwait,boolean *Uwait,
				int *DmaxDis, int *UmaxDis,int maxNodes)
{	
	if(*Dwait)
		return 1;

	unsigned int currNode,twin,tempNode;
	CTGinHEAP *ctgInHeap;
	int indexInArray;
	CONNECT *us_cnt,*ds_cnt;
	int dis0,dis;
	boolean ret,isEmpty;

	while((indexInArray=removeNextNodeFromHeap(dheap))!=0){
		ctgInHeap = &ctg4heapArray[indexInArray];
		currNode = ctgInHeap->ctgID;
		dis0 = ctgInHeap->dis;

		isEmpty = IsHeapEmpty(dheap);

		twin = getTwinCtg(currNode);
		us_cnt = ctgInHeap->us_shut4dheap? NULL:contig_array[twin].downwardConnect;
		while(us_cnt){
			if(us_cnt->deleted||us_cnt->mask||
				contig_array[getTwinCtg(us_cnt->contigID)].inSubGraph){
				us_cnt = us_cnt->next;
				continue;
			}	

			tempNode = getTwinCtg(us_cnt->contigID);
			if(contig_array[tempNode].inSubGraph){
				us_cnt = us_cnt->next;
				continue;
			}
			dis = dis0 - us_cnt->gapLen - (int)contig_array[twin].length;

			ret = dispatch1node(dis,tempNode,maxNodes,dheap,uheap,DmaxDis,UmaxDis);
			if(ret==0)
				return 0;	
			else if(ret<0)
				*Uwait = 0;
			
			us_cnt = us_cnt->next;
		}
		
		if(nodeCounter>1&&isEmpty){
			*Dwait = canDheapWait(currNode,dis0,*DmaxDis);	
			if(*Dwait){
				isEmpty = IsHeapEmpty(dheap);
				insertNodeIntoHeap(dheap,dis0,indexInArray);
				ctg4heapArray[indexInArray].us_shut4dheap = 1;
				if(isEmpty)
					return 1;
				else 
					continue;
			}
		}
		ds_cnt = ctgInHeap->ds_shut4dheap? NULL:contig_array[currNode].downwardConnect;
		while(ds_cnt){
			if(ds_cnt->deleted||ds_cnt->mask||contig_array[ds_cnt->contigID].inSubGraph){
				ds_cnt = ds_cnt->next;
				continue;
			}	
			tempNode = ds_cnt->contigID;
			dis = dis0 + ds_cnt->gapLen + (int)contig_array[tempNode].length;
			ret = dispatch1node(dis,tempNode,maxNodes,dheap,uheap,DmaxDis,UmaxDis);
			if(ret==0)
				return 0;	
			else if(ret<0)
				*Uwait = 0;
		}  // for each downstream connections
	}  // for each node comes off the heap
	
	*Dwait = 1;
	return 1;
}

static boolean canUheapWait(unsigned int currNode,int dis, int UmaxDis)
{
	int temp_len = contig_array[currNode].length;
	if(-dis+temp_len<UmaxDis)
		return 0;
	else
		return 1;
}

static boolean workOnUheap(FibHeap *dheap,FibHeap *uheap,boolean *Dwait, boolean *Uwait,
				int *DmaxDis, int *UmaxDis,int maxNodes)
{
	if(*Uwait)
		return 1;
	unsigned int currNode,twin,tempNode;
	CTGinHEAP *ctgInHeap;
	int indexInArray;
	CONNECT *us_cnt,*ds_cnt;
	int dis0,dis;
	boolean ret,isEmpty;

	while((indexInArray=removeNextNodeFromHeap(uheap))!=0){
		ctgInHeap = &ctg4heapArray[indexInArray];
		currNode = ctgInHeap->ctgID;
		dis0 = ctgInHeap->dis;

		isEmpty = IsHeapEmpty(uheap);
		ds_cnt = ctgInHeap->ds_shut4uheap? NULL:contig_array[currNode].downwardConnect;
		while(ds_cnt){
			if(ds_cnt->deleted||ds_cnt->mask||contig_array[ds_cnt->contigID].inSubGraph){
				ds_cnt = ds_cnt->next;
				continue;
			}	
			tempNode = ds_cnt->contigID;
			dis = dis0 + ds_cnt->gapLen + contig_array[tempNode].length;
			ret = dispatch1node(dis,tempNode,maxNodes,dheap,uheap,DmaxDis,UmaxDis);
			if(ret==0)
				return 0;	
			else if(ret>0)
				*Dwait = 0;
			
		}  // for each downstream connections

		if(nodeCounter>1&&isEmpty){
			*Uwait = canUheapWait(currNode,dis0,*UmaxDis);	
			if(*Uwait){
				isEmpty = IsHeapEmpty(uheap);
				insertNodeIntoHeap(uheap,dis0,indexInArray);
				ctg4heapArray[indexInArray].ds_shut4uheap = 1;
				if(isEmpty)
					return 1;
				else
					continue;
			}
		}

		twin = getTwinCtg(currNode);
		us_cnt = ctgInHeap->us_shut4uheap? NULL:contig_array[twin].downwardConnect;
		while(us_cnt){
			if(us_cnt->deleted||us_cnt->mask||
				contig_array[getTwinCtg(us_cnt->contigID)].inSubGraph){
				us_cnt = us_cnt->next;
				continue;
			}	

			tempNode = getTwinCtg(us_cnt->contigID);
			if(contig_array[tempNode].inSubGraph){
				us_cnt = us_cnt->next;
				continue;
			}
			dis = dis0 - us_cnt->gapLen - contig_array[twin].length;

			ret = dispatch1node(dis,tempNode,maxNodes,dheap,uheap,DmaxDis,UmaxDis);
			if(ret==0)
				return 0;	
			else if(ret>0)
				*Dwait = 1;
			
			us_cnt = us_cnt->next;
		}
		
	}  // for each node comes off the heap
	
	*Uwait = 1;
	return 1;
}

static boolean pickUpGeneralSubgraph(unsigned int node1,int maxNodes)
{
	FibHeap *Uheap = newFibHeap();  // heap for upstream contigs to node1
	FibHeap *Dheap = newFibHeap();
	int UmaxDis;   // max distance upstream to node1 
	int DmaxDis;
	boolean Uwait;  // wait signal for Uheap
	boolean Dwait;
	int dis;
	boolean ret;
	
	//initiate: node1 is put to array once, and to both Dheap and Uheap
	dis = 0;
	nodeCounter = 1;
	putNodeInArray(node1,maxNodes,dis);
	insertNodeIntoHeap(Dheap,dis,nodeCounter);
	ctg4heapArray[nodeCounter].us_shut4dheap = 1;
	Dwait = 0;
	DmaxDis = 0;
	
	insertNodeIntoHeap(Uheap,dis,nodeCounter);
	ctg4heapArray[nodeCounter].ds_shut4uheap = 1;
	Uwait = 1;
	UmaxDis = contig_array[node1].length;

	while(1){
		ret = workOnDheap(Dheap,Uheap,&Dwait,&Uwait,&DmaxDis,&UmaxDis,maxNodes);
		if(!ret){
 			setInGraph(0);
			destroyHeap(Dheap);
			destroyHeap(Uheap);
			return 0;
		}
		ret = workOnUheap(Dheap,Uheap,&Dwait,&Uwait,&DmaxDis,&UmaxDis,maxNodes);
		if(!ret){
 			setInGraph(0);
			destroyHeap(Dheap);
			destroyHeap(Uheap);
			return 0;
		}
		if(Uwait&&Dwait){
			destroyHeap(Dheap);
			destroyHeap(Uheap);
			return 1;
		}
	}
	
}

static int cmp_ctg(const void *a,const void *b)
{
	CTGinHEAP *A,*B;
	A = (CTGinHEAP *)a;
	B = (CTGinHEAP *)b;

	if(A->dis>B->dis)
		return 1;
	else if(A->dis==B->dis)
		return 0;
	else
		return -1;
}

static boolean checkEligible()
{
	unsigned int firstNode = ctg4heapArray[1].ctgID;
	unsigned int twin;
	int i;
	boolean flag = 0;
	
	//check if the first node has incoming link from twin of any node in subgraph
	// or it has multi outgoing links bound to incoming links
	twin = getTwinCtg(firstNode);
	CONNECT *ite_cnt = contig_array[twin].downwardConnect;
	while(ite_cnt){
		if(ite_cnt->deleted||ite_cnt->mask){
			ite_cnt = ite_cnt->next;
			continue;
		}
		if(contig_array[ite_cnt->contigID].inSubGraph){
/*
			if(firstNode==3693)
				printf("eligible link %d -> %d\n",twin,ite_cnt->contigID);
*/
			return 0;
		}
		if(ite_cnt->prevInScaf){
			if(flag)  
				return 0;
			flag = 1;
		}
		ite_cnt = ite_cnt->next;
	}

	//check if the last node has outgoing link to twin of any node in subgraph
	// or it has multi outgoing links bound to incoming links
	unsigned int lastNode = ctg4heapArray[nodeCounter].ctgID;
	ite_cnt = contig_array[lastNode].downwardConnect;
	flag = 0;
	while(ite_cnt){
		if(ite_cnt->deleted||ite_cnt->mask){
			ite_cnt = ite_cnt->next;
			continue;
		}
		twin = getTwinCtg(ite_cnt->contigID);
		if(contig_array[twin].inSubGraph){
/*
			if(firstNode==3693)
				printf("eligible link %d -> %d\n",lastNode,ite_cnt->contigID);
*/
			return 0;
		}
		if(ite_cnt->prevInScaf){
			if(flag)  
				return 0;
			flag = 1;
		}
		ite_cnt = ite_cnt->next;
	}
	//check if any node has outgoing link to node outside the subgraph
	for(i=1;i<nodeCounter;i++){
		ite_cnt = contig_array[ctg4heapArray[i].ctgID].downwardConnect;
		while(ite_cnt){
			if(ite_cnt->deleted||ite_cnt->mask){
				ite_cnt = ite_cnt->next;
				continue;
			}
			if(!contig_array[ite_cnt->contigID].inSubGraph){
			/*
				printf("eligible check: ctg %d links to ctg %d\n",
					ctg4heapArray[i].ctgID,ite_cnt->contigID);
			*/
				return 0;
			}
			ite_cnt = ite_cnt->next;
		}
	}
	//check if any node has incoming link from node outside the subgraph
	for(i=2;i<=nodeCounter;i++){
		twin = getTwinCtg(ctg4heapArray[i].ctgID);
		ite_cnt = contig_array[twin].downwardConnect;
		while(ite_cnt){
			if(ite_cnt->deleted||ite_cnt->mask){
				ite_cnt = ite_cnt->next;
				continue;
			}
			if(!contig_array[getTwinCtg(ite_cnt->contigID)].inSubGraph){
			/*
				printf("eligible check: ctg %d links to ctg %d\n",
					ctg4heapArray[i].ctgID,ite_cnt->contigID);
			*/
				return 0;
			}
			ite_cnt = ite_cnt->next;
		}
	}

	return 1;
}

//put nodes in sub-graph in a line
static void arrangeNodes_general()
{
	int i,gap;
	CONNECT *ite_cnt,*temp_cnt,*bal_cnt,*prev_cnt,*next_cnt;
	unsigned int node1,node2;
	unsigned int bal_nd1,bal_nd2;
	//delete original connections
	for(i=1;i<=nodeCounter;i++){
		node1 = ctg4heapArray[i].ctgID;
		ite_cnt = contig_array[node1].downwardConnect;
		while(ite_cnt){
			if(ite_cnt->mask||ite_cnt->deleted||!contig_array[ite_cnt->contigID].inSubGraph){
				ite_cnt = ite_cnt->next;
				continue;
			}
			ite_cnt->deleted = 1;
			setNextInScaf(ite_cnt,NULL);
			setPrevInScaf(ite_cnt,0);
			ite_cnt = ite_cnt->next;
		}

		bal_nd1 = getTwinCtg(node1);
		ite_cnt = contig_array[bal_nd1].downwardConnect;
		while(ite_cnt){
			if(ite_cnt->mask||ite_cnt->deleted||!contig_array[getTwinCtg(ite_cnt->contigID)].inSubGraph){
				ite_cnt = ite_cnt->next;
				continue;
			}
			ite_cnt->deleted = 1;
			setNextInScaf(ite_cnt,NULL);
			setPrevInScaf(ite_cnt,0);
			ite_cnt = ite_cnt->next;
		}
	}
	//create new connections
	prev_cnt = next_cnt = NULL;
	for(i=1;i<nodeCounter;i++){
		node1 = ctg4heapArray[i].ctgID;
		node2 = ctg4heapArray[i+1].ctgID;
		bal_nd1 = getTwinCtg(node1);
		bal_nd2 = getTwinCtg(node2);
		gap = ctg4heapArray[i+1].dis - ctg4heapArray[i].dis
				- contig_array[node2].length;
		temp_cnt = getCntBetween(node1,node2);
		if(temp_cnt){
			temp_cnt->deleted = 0;
			temp_cnt->mask = 0;
			//temp_cnt->gapLen = gap;
			bal_cnt = getCntBetween(bal_nd2,bal_nd1);
			bal_cnt->deleted = 0;
			bal_cnt->mask = 0;
			//bal_cnt->gapLen = gap;
		}
		else{
			temp_cnt = allocateCN(node2,gap);
			if(cntLookupTable)
				putCnt2LookupTable(node1,temp_cnt);
			temp_cnt->next = contig_array[node1].downwardConnect;
			contig_array[node1].downwardConnect = temp_cnt;
			bal_cnt = allocateCN(bal_nd1,gap);
			if(cntLookupTable)
				putCnt2LookupTable(bal_nd2,bal_cnt);
			bal_cnt->next = contig_array[bal_nd2].downwardConnect;
			contig_array[bal_nd2].downwardConnect = bal_cnt;
		}
		if(prev_cnt){
			setNextInScaf(prev_cnt,temp_cnt);
			setPrevInScaf(temp_cnt,1);
		}
		if(next_cnt){
			setNextInScaf(bal_cnt,next_cnt);
			setPrevInScaf(next_cnt,1);
		}
		prev_cnt = temp_cnt;
		next_cnt = bal_cnt;
	}

	//re-binding connection at both ends
	bal_nd2 = getTwinCtg(ctg4heapArray[1].ctgID);
	ite_cnt = contig_array[bal_nd2].downwardConnect;
	while(ite_cnt){
		if(ite_cnt->deleted||ite_cnt->mask){
			ite_cnt = ite_cnt->next;
			continue;
		}
		if(ite_cnt->prevInScaf)
			break;
		ite_cnt = ite_cnt->next;
	}
	if(ite_cnt){
		bal_nd1 = ite_cnt->contigID;
		node1 = getTwinCtg(bal_nd1);
		node2 = ctg4heapArray[1].ctgID;
		temp_cnt = checkConnect(node1,node2);
		bal_cnt = ite_cnt;
		next_cnt = checkConnect(ctg4heapArray[1].ctgID,ctg4heapArray[2].ctgID);
		prev_cnt = checkConnect(getTwinCtg(ctg4heapArray[2].ctgID), getTwinCtg(ctg4heapArray[1].ctgID));
		if(temp_cnt){
			setNextInScaf(temp_cnt,next_cnt);
			setPrevInScaf(temp_cnt->nextInScaf,0);
			setPrevInScaf(next_cnt,1);
			setNextInScaf(prev_cnt,bal_cnt);
		}
	}

	node1 = ctg4heapArray[nodeCounter].ctgID;
	ite_cnt = contig_array[node1].downwardConnect;
	while(ite_cnt){
		if(ite_cnt->deleted||ite_cnt->mask){
			ite_cnt = ite_cnt->next;
			continue;
		}
		if(ite_cnt->prevInScaf)
			break;
		ite_cnt = ite_cnt->next;
	}
	if(ite_cnt){
		node2 = ite_cnt->contigID;
		bal_nd1 = getTwinCtg(node1);
		bal_nd2 = getTwinCtg(node2);
		temp_cnt = ite_cnt;
		bal_cnt = checkConnect(bal_nd2,bal_nd1);
		next_cnt = checkConnect(getTwinCtg(ctg4heapArray[nodeCounter].ctgID),
					getTwinCtg(ctg4heapArray[nodeCounter-1].ctgID));
		prev_cnt = checkConnect(ctg4heapArray[nodeCounter-1].ctgID,ctg4heapArray[nodeCounter].ctgID);
		setNextInScaf(prev_cnt,temp_cnt);
		setNextInScaf(bal_cnt,next_cnt);
		setPrevInScaf(next_cnt,1);
	}
}
//check if contigs next to each other have reasonable overlap
boolean checkOverlapInBetween_general(double tolerance)
{
	int i,gap;
	unsigned int node;
	int lenSum,lenOlp;
	lenSum = lenOlp = 0;
	for(i=1;i<=nodeCounter;i++){
		node = ctg4heapArray[i].ctgID;
		lenSum += contig_array[node].length;
	}
	if(lenSum<1)
		return 1;
	for(i=1;i<nodeCounter;i++){
		gap = ctg4heapArray[i+1].dis - ctg4heapArray[i].dis
				- contig_array[ctg4heapArray[i+1].ctgID].length;
		if(-gap>0)
			lenOlp += -gap;
		//if(-gap>ins_size_var)

	}
	double olp_pect=(double)lenOlp/lenSum;
	fprintf(stderr,"[%s]existing with olp_pect %.3f.\n",__FUNCTION__,olp_pect);
	if(olp_pect>tolerance){
		return 0;
	}
	return 1;
}

//check if there's any connect indicates the opposite order between nodes in sub-graph
static boolean checkConflictCnt_general(double tolerance)
{
	int i,j;
	int supportCounter=0;
	int objectCounter=0;
	CONNECT *cnt;
	for(i=1;i<nodeCounter;i++){
		for(j=i+1;j<=nodeCounter;j++){	
			//cnt=getCntBetween(nodesInSubInOrder[j],nodesInSubInOrder[i]);
			cnt = checkConnect(ctg4heapArray[i].ctgID,ctg4heapArray[j].ctgID);
			if(cnt)
				supportCounter += cnt->weight;
			cnt = checkConnect(ctg4heapArray[j].ctgID,ctg4heapArray[i].ctgID);
			if(cnt)
				objectCounter += cnt->weight;
				//return 1;
		}
	}
	if(supportCounter<1)
		return 1;
	if((double)objectCounter/supportCounter<tolerance)
		return 0;
	return 1;
}
// turn sub-graph to linear struct
static void general_linearization(boolean strict)
{
	unsigned int i;
	int subCounter=0;
	int out_num;
	boolean flag;
	int conflCounter=0,overlapCounter=0,eligibleCounter=0;
	double overlapTolerance,conflTolerance;
	
	for(i=num_ctg;i>0;i--){
		if(contig_array[i].mask)
			continue;
		out_num = validConnect(i,NULL);

		if(out_num<2)
			continue;
		
		//flag = pickSubGraph(i,strict);
		flag = pickUpGeneralSubgraph(i,MaxNodeInSub);
		if(!flag)
			continue;
		subCounter++;
		qsort(&ctg4heapArray[1],nodeCounter,sizeof(CTGinHEAP),cmp_ctg);
		flag = checkEligible();
		if(!flag){
			eligibleCounter++;
 			setInGraph(0);
			continue;
		}
		if(strict){
			overlapTolerance = OverlapPercent;
			conflTolerance = ConflPercent;
		}else{
			overlapTolerance = 2*OverlapPercent;
			conflTolerance = 2*ConflPercent;
		}
		flag = checkOverlapInBetween_general(overlapTolerance);
		if(!flag){
			overlapCounter++;
 			setInGraph(0);
			continue;
		}
		flag = checkConflictCnt_general(conflTolerance);
		if(flag){
			conflCounter++;
 			setInGraph(0);
			continue;
		}
		arrangeNodes_general();
 		setInGraph(0);
	}
	fprintf(stdout,"[%s]Picked  %d subgraphs,%d have conflicting connections,%d have significant overlapping, %d eligible\n",
			__FUNCTION__,subCounter,conflCounter,overlapCounter,eligibleCounter);

}

/****       the fowllowing codes for detecting and break down scaffold at weak point  **********/
// mark connections in scaffolds made by small pe 
static void smallScaf()
{
	unsigned int i,ctg,bal_ctg,prevCtg;
	int counter=0;
	CONNECT *bindCnt,*cnt;

	for(i=1;i<=num_ctg;i++)
		contig_array[i].flag = 0;
	for(i=1;i<=num_ctg;i++){
		if(contig_array[i].flag||contig_array[i].mask||!contig_array[i].downwardConnect)
			continue;
		bindCnt = getBindCnt(i);
		if(!bindCnt)
			continue;
		counter++;
		
		contig_array[i].flag = 1;
		contig_array[getTwinCtg(i)].flag = 1;
		prevCtg = getTwinCtg(i);
		while(bindCnt){
			ctg = bindCnt->contigID;
			bal_ctg = getTwinCtg(ctg);
			bindCnt->bySmall = 1;
			cnt = getCntBetween(bal_ctg,prevCtg);
			if(cnt)
				cnt->bySmall = 1;
			
			contig_array[ctg].flag = 1;
			contig_array[bal_ctg].flag = 1;
			prevCtg = bal_ctg;
			bindCnt = bindCnt->nextInScaf;
		}

		ctg = getTwinCtg(i);
		bindCnt = getBindCnt(ctg);
		prevCtg = i;
		while(bindCnt){
			ctg = bindCnt->contigID;
			bal_ctg = getTwinCtg(ctg);
			bindCnt->bySmall = 1;
			cnt = getCntBetween(bal_ctg,prevCtg);
			if(cnt)
				cnt->bySmall = 1;
			
			contig_array[ctg].flag = 1;
			contig_array[bal_ctg].flag = 1;
			prevCtg = bal_ctg;
			bindCnt = bindCnt->nextInScaf;
		}
	}
	//printf("Report from smallScaf: %d scaffolds by smallPE\n",counter);
}

static boolean putItem2Sarray(unsigned int scaf,int wt,DARRAY *SCAF,DARRAY *WT,int counter)
{
	int i;
	unsigned int *scafP,*wtP;
	for(i=0;i<counter;i++){
		scafP = (unsigned int *)darrayGet(SCAF,i);
		if((*scafP)==scaf){
			wtP = (unsigned int *)darrayGet(WT,i);
			*wtP = (*wtP + wt);
			return 0;
		}
	}
	scafP = (unsigned int *)darrayPut(SCAF,counter);
	wtP = (unsigned int *)darrayPut(WT,counter);
	*scafP = scaf;
	*wtP = wt;
	return 1;
}

static int getDSLink2Scaf(STACK *scafStack,DARRAY *SCAF,DARRAY *WT)
{
	CONNECT *ite_cnt;
	unsigned int ctg,targetCtg,*pt;
	int counter=0;
	boolean inc;
	
	stackRecover(scafStack);
	while((pt=(unsigned int*)stackPop(scafStack))!=NULL){
		ctg = *pt;
		if(contig_array[ctg].mask||!contig_array[ctg].downwardConnect)
			continue;
		ite_cnt = contig_array[ctg].downwardConnect;
		while(ite_cnt){
			if(ite_cnt->deleted||ite_cnt->mask||ite_cnt->singleInScaf
				||ite_cnt->nextInScaf||ite_cnt->prevInScaf||ite_cnt->inherit){
				ite_cnt = ite_cnt->next;	
				continue;
			}
			targetCtg = ite_cnt->contigID;
			if(contig_array[ctg].from_vt==contig_array[targetCtg].from_vt){  // on the same scaff
				ite_cnt = ite_cnt->next;	
				continue;
			}
			inc = putItem2Sarray(contig_array[targetCtg].from_vt,ite_cnt->weight,SCAF,WT,counter);
			if(inc)
				counter++;
			ite_cnt = ite_cnt->next;	
		}
	}
	return counter;

}

static int getScaffold(unsigned int start, STACK *scafStack)
{
	int len = contig_array[start].length;
	unsigned int *pt,ctg;

	emptyStack(scafStack);
	pt = (unsigned int*)stackPush(scafStack);
	*pt = start;
	CONNECT *bindCnt = getBindCnt(start);
	while(bindCnt){
		ctg = bindCnt->contigID;
		pt = (unsigned int*)stackPush(scafStack);
		*pt = ctg;
		len += contig_array[ctg].length;
		bindCnt = bindCnt->nextInScaf;
	}
	stackBackup(scafStack);
	return len;
}

static boolean isLinkReliable(DARRAY *WT,int count)
{
	int i;
	for(i=0;i<count;i++)
		if(*(int *)darrayGet(WT,i)>=weakPE)
			return 1;

	return 0;
}

static int getWtFromSarray(DARRAY *SCAF,DARRAY *WT,int count,unsigned int scaf)
{
	int i;
	for(i=0;i<count;i++)
		if(*(unsigned int *)darrayGet(SCAF,i)==scaf)
			return *(int *)darrayGet(WT,i);

	return 0;
}

static void switch2twin(STACK *scafStack)
{
	unsigned int *pt;
	stackRecover(scafStack);
	while((pt=(unsigned int*)stackPop(scafStack))!=NULL)
		*pt = getTwinCtg(*pt);
}
/*
           ------>
 scaf1 --- --- -- -- --- 
                         scaf2 -- --- --- --
                                   ---->
*/
static boolean checkScafConsist(STACK *scafStack1,STACK *scafStack2)
{
	DARRAY *downwardTo1 = (DARRAY *)createDarray(1000,sizeof(unsigned int));// scaf links to those scaffolds
	DARRAY *downwardTo2 = (DARRAY *)createDarray(1000,sizeof(unsigned int));
	DARRAY *downwardWt1 = (DARRAY *)createDarray(1000,sizeof(unsigned int));// scaf links to scaffolds with those wt
	DARRAY *downwardWt2 = (DARRAY *)createDarray(1000,sizeof(unsigned int));

	int linkCount1 = getDSLink2Scaf(scafStack1,downwardTo1,downwardWt1);
	int linkCount2 = getDSLink2Scaf(scafStack2,downwardTo2,downwardWt2);
	if(!linkCount1||!linkCount2){
		freeDarray(downwardTo1);
		freeDarray(downwardTo2);
		freeDarray(downwardWt1);
		freeDarray(downwardWt2);
		return 1; 
	}
	boolean flag1 = isLinkReliable(downwardWt1,linkCount1);
	boolean flag2 = isLinkReliable(downwardWt2,linkCount2);
	if(!flag1||!flag2){
		freeDarray(downwardTo1);
		freeDarray(downwardTo2);
		freeDarray(downwardWt1);
		freeDarray(downwardWt2);
		return 1; 
	}

	unsigned int scaf;
	int i,wt1,wt2,ret=1;

	for(i=0;i<linkCount1;i++){
		wt1 = *(int *)darrayGet(downwardWt1,i);
		if(wt1<weakPE)
			continue;
		scaf = *(unsigned int *)darrayGet(downwardTo1,i);
		wt2 = getWtFromSarray(downwardTo2,downwardWt2,linkCount2,scaf);
		if(wt2<1){
			//fprintf(stderr,"Inconsistant link to %d\n",scaf);
			ret = 0;
			break;
		}
	}

	freeDarray(downwardTo1);
	freeDarray(downwardTo2);
	freeDarray(downwardWt1);
	freeDarray(downwardWt2);
	return ret;
}

static void setBreakPoints(DARRAY *ctgArray,int count, int weakest,
				int *start,int *finish) 
{
	int index=weakest-1;
	unsigned int thisCtg;
	unsigned int nextCtg = *(unsigned int *)darrayGet(ctgArray,weakest);
	CONNECT *cnt;
	*start = weakest;
	while(index>=0){
		thisCtg = *(unsigned int *)darrayGet(ctgArray,index);
		cnt = getCntBetween(thisCtg,nextCtg);
		if(cnt->maxGap>2)
			break;
		else
			*start = index;	
		nextCtg = thisCtg;
		index--;
	}
	unsigned int prevCtg = *(unsigned int *)darrayGet(ctgArray,weakest+1);
	*finish = weakest+1;
	index = weakest+2;
	while(index<count){
		thisCtg = *(unsigned int *)darrayGet(ctgArray,index);
		cnt = getCntBetween(prevCtg,thisCtg);
		if(cnt->maxGap>2)
			break;
		else
			*finish = index;	
		prevCtg = thisCtg;
		index++;
	}

}

static void changeScafEnd(STACK *scafStack,unsigned int end)	
{
	
	unsigned int ctg,*pt;
	unsigned int start=getTwinCtg(end);
	stackRecover(scafStack);
	while((pt=(unsigned int*)stackPop(scafStack))!=NULL){
		ctg = *pt;
		contig_array[ctg].to_vt = end;
		contig_array[getTwinCtg(ctg)].from_vt = start;
	}
}

static void changeScafBegin(STACK *scafStack,unsigned int start)	
{
	
	unsigned int ctg,*pt;
	unsigned int end=getTwinCtg(start);
	stackRecover(scafStack);
	while((pt=(unsigned int*)stackPop(scafStack))!=NULL){
		ctg = *pt;
		contig_array[ctg].from_vt = start;
		contig_array[getTwinCtg(ctg)].to_vt = end;
	}
}
// break down scaffolds poorly supported by longer PE
static void detectBreakScaf()
{
	fprintf(stderr,"[%s]entering this function.\n",__FUNCTION__);
	unsigned int i,avgPE,scafLen,len,ctg,bal_ctg,prevCtg,thisCtg;
	long long peCounter,linkCounter;
	int num3,num5,weakPoint,tempCounter,j,t,counter=0;
	CONNECT *bindCnt,*cnt,*weakCnt;

	STACK *scafStack1 = (STACK *)createStack(1000,sizeof(unsigned int));
	STACK *scafStack2 = (STACK *)createStack(1000,sizeof(unsigned int));

	for(i=1;i<=num_ctg;i++)
		contig_array[i].flag = 0;
	for(i=1;i<=num_ctg;i++){
		if(contig_array[i].flag||contig_array[i].mask||!contig_array[i].downwardConnect)
			continue;
		bindCnt = getBindCnt(i);
		if(!bindCnt)
			continue;
		//first scan get the average coverage by longer pe
		num5 = num3 = peCounter = linkCounter = 0;
		scafLen = contig_array[i].length;
		ctg = i;
		*(unsigned int *)darrayPut(scaf5,num5++) = i;
		contig_array[i].flag = 1;
		contig_array[getTwinCtg(i)].flag = 1;
		while(bindCnt){
			if(!bindCnt->bySmall)
				break;
			linkCounter++;
			peCounter += bindCnt->maxGap;
			ctg = bindCnt->contigID;
			scafLen += contig_array[ctg].length;
			*(unsigned int *)darrayPut(scaf5,num5++) = ctg;
			bal_ctg = getTwinCtg(ctg);
			contig_array[ctg].flag = 1;
			contig_array[bal_ctg].flag = 1;
			bindCnt = bindCnt->nextInScaf;
		}

		ctg = getTwinCtg(i);
		bindCnt = getBindCnt(ctg);
		while(bindCnt){
			if(!bindCnt->bySmall)
				break;
			linkCounter++;
			peCounter += bindCnt->maxGap;
			ctg = bindCnt->contigID;
			scafLen += contig_array[ctg].length;
			bal_ctg = getTwinCtg(ctg);
			contig_array[ctg].flag = 1;
			contig_array[bal_ctg].flag = 1;
			*(unsigned int *)darrayPut(scaf3,num3++) = bal_ctg;
			bindCnt = bindCnt->nextInScaf;
		}
		if(linkCounter<1||scafLen<5000)
			continue;
	
		avgPE = peCounter/linkCounter;
		
		if(avgPE<10)
			continue;

		tempCounter = 0;
		for(j=num3-1;j>=0;j--)
			*(unsigned int *)darrayPut(tempArray,tempCounter++) = 	
						*(unsigned int *)darrayGet(scaf3,j);
		
		for(j=0;j<num5;j++)
			*(unsigned int *)darrayPut(tempArray,tempCounter++) = 
					*(unsigned int *)darrayGet(scaf5,j);
		
		prevCtg = *(unsigned int *)darrayGet(tempArray,0);
		weakCnt = NULL;
		weakPoint = 0;
		len = contig_array[prevCtg].length;
		for(t=1;t<tempCounter;t++){
			thisCtg = *(unsigned int *)darrayGet(tempArray,t);		
			if(len<2000){
				len += contig_array[thisCtg].length;
				prevCtg = thisCtg;
				continue;
			}else if(len>scafLen-2000)
				break;
			len += contig_array[thisCtg].length;
			if(contig_array[prevCtg].from_vt!=contig_array[thisCtg].from_vt||
				contig_array[prevCtg].indexInScaf>contig_array[thisCtg].indexInScaf){
				prevCtg = thisCtg;
				continue;
			}
			cnt = getCntBetween(prevCtg,thisCtg);
			if(!weakCnt||weakCnt->maxGap>cnt->maxGap){
				weakCnt = cnt;
				weakPoint = t;
			}
			prevCtg = thisCtg;
		}
		if(!weakCnt||(weakCnt->maxGap>2&&weakCnt->maxGap>avgPE/5))
			continue;
		prevCtg = *(unsigned int *)darrayGet(tempArray,weakPoint-1);
		thisCtg = *(unsigned int *)darrayGet(tempArray,weakPoint);
		if(contig_array[prevCtg].from_vt!=contig_array[thisCtg].from_vt||
			contig_array[prevCtg].indexInScaf>contig_array[thisCtg].indexInScaf){
			printf("contig %d and %d not on the same scaff\n",prevCtg,thisCtg);
			continue;
		}
		setConnectWP(prevCtg,thisCtg,1);
		/*
		fprintf(stderr,"scaffold len %d, avg long pe cov %d (%ld/%ld)\n",
					scafLen,avgPE,peCounter,linkCounter);
		fprintf(stderr,"Weak connect (%d) between %d(%dth of %d) and %d\n"
					,weakCnt->maxGap,prevCtg,weakPoint-1,tempCounter,thisCtg);
		*/
		// set start and end to break down the scaffold
		int index1,index2;
		setBreakPoints(tempArray,tempCounter,weakPoint-1,&index1,&index2); 
		//fprintf(stderr,"break %d ->...-> %d\n",index1,index2);
		unsigned int start = *(unsigned int*)darrayGet(tempArray,index1);
		unsigned int finish = *(unsigned int*)darrayGet(tempArray,index2);
		int len1 = getScaffold(getTwinCtg(start), scafStack1);
		int len2 = getScaffold(finish, scafStack2);
		if(len1<2000||len2<2000)
			continue;
		switch2twin(scafStack1);
		int flag1 = checkScafConsist(scafStack1,scafStack2);

		switch2twin(scafStack1);
		switch2twin(scafStack2);
		int flag2 = checkScafConsist(scafStack2,scafStack1);
		if(!flag1||!flag2){
			changeScafBegin(scafStack1,getTwinCtg(start));	
			changeScafEnd(scafStack2,getTwinCtg(finish));	
			//unbind links
			unsigned int nextCtg = *(unsigned int *)darrayGet(tempArray,index1+1);
			thisCtg = *(unsigned int *)darrayGet(tempArray,index1);
			cnt=getCntBetween(getTwinCtg(nextCtg),getTwinCtg(thisCtg));
			if(cnt->nextInScaf){
				prevCtg = getTwinCtg(cnt->nextInScaf->contigID);
				cnt->nextInScaf->prevInScaf = 0;
				cnt = getCntBetween(prevCtg,thisCtg);
				cnt->nextInScaf = NULL;
			}
			prevCtg = *(unsigned int *)darrayGet(tempArray,index2-1);
			thisCtg = *(unsigned int *)darrayGet(tempArray,index2);
			cnt = getCntBetween(prevCtg,thisCtg);
			if(cnt->nextInScaf){
				nextCtg = cnt->nextInScaf->contigID;
				cnt->nextInScaf->prevInScaf= 0;
				cnt = getCntBetween(getTwinCtg(nextCtg),getTwinCtg(thisCtg));
				cnt->nextInScaf = NULL;
			}
			prevCtg = *(unsigned int *)darrayGet(tempArray,index1);
			for(t=index1+1;t<=index2;t++){
				thisCtg = *(unsigned int *)darrayGet(tempArray,t);		
				cnt = getCntBetween(prevCtg,thisCtg);
				cnt->mask = 1;
				cnt->nextInScaf=NULL;
				cnt->prevInScaf = 0;
				cnt = getCntBetween(getTwinCtg(thisCtg),getTwinCtg(prevCtg));
				cnt->mask = 1;
				cnt->nextInScaf=NULL;
				cnt->prevInScaf = 0;
				/*
				fprintf(stderr,"(%d %d)/(%d %d) ",
					prevCtg,thisCtg,getTwinCtg(thisCtg),getTwinCtg(prevCtg));
				*/
				prevCtg = thisCtg;
			}
			//fprintf(stderr,": BREAKING\n");
			counter++;
		}
	}

	freeStack(scafStack1);
	freeStack(scafStack2);
	fprintf(stderr,"[%s]existing this function.\n",__FUNCTION__);
	//printf("Report from checkScaf: %d scaffold segments broken\n",counter);
}

static boolean checkSimple(DARRAY *ctgArray,int count)
{
	int i;
	unsigned int ctg;
	for(i=0;i<count;i++){
		ctg = *(unsigned int *)darrayGet(ctgArray,i);
		contig_array[ctg].flag = 0;
		contig_array[getTwinCtg(ctg)].flag = 0;
	}
	for(i=0;i<count;i++){
		ctg = *(unsigned int *)darrayGet(ctgArray,i);
		if(contig_array[ctg].flag)
			return 0;
		contig_array[ctg].flag = 1;
		contig_array[getTwinCtg(ctg)].flag = 1;
	}
	return 1;
		
}

static void checkCircle()
{
	unsigned int i,ctg;
	CONNECT *cn_temp1;
	int counter=0;

	for(i=1;i<=num_ctg;i++){
		cn_temp1 = contig_array[i].downwardConnect;
		while(cn_temp1){
			if(cn_temp1->weak||cn_temp1->deleted){
				cn_temp1 = cn_temp1->next;
				continue;
			}
			ctg = cn_temp1->contigID;
			if(checkConnect(ctg,i)){
				counter++;
				maskContig(i,1);
				maskContig(ctg,1);
			}
			cn_temp1 = cn_temp1->next;
		}

	}
	//printf("%d circles removed \n",counter);
}

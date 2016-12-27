#include "stdinc.h" 
#include "newhash.h"
#include "extfunc.h" 
#include "extvab.h" 

static int Ncounter;
static int allGaps;

// for multi threads
static int scafBufSize = 100;
static STACK **ctgStackBuffer;
static int scafCounter;
static int scafInBuf;

static void convertIndex()
{
	int *length_array = (int *)ckalloc((num_ctg+1)*sizeof(int));
	unsigned int i;
	for(i=1;i<=num_ctg;i++)
		length_array[i] = 0;

	for(i=1;i<=num_ctg;i++){
		if(index_array[i]>0)
			length_array[index_array[i]] = i;
	}
	for(i=1;i<=num_ctg;i++)
		index_array[i] = length_array[i];  //contig i with new index: index_array[i]
	free((void *)length_array);

}

static void reverseStack(STACK *dStack,STACK *sStack)
{
	CTGinSCAF *actg,*ctgPt;
	emptyStack(dStack);

	while((actg=(CTGinSCAF *)stackPop(sStack))!=NULL){
		ctgPt = (CTGinSCAF *)stackPush(dStack);
		ctgPt->ctgID = actg->ctgID;
		ctgPt->start = actg->start;
		ctgPt->end = actg->end;
	}
	stackBackup(dStack);
}

static void initStackBuf(STACK **ctgStackBuffer,int scafBufSize)
{
	int i;
	for(i=0;i<scafBufSize;i++)
		ctgStackBuffer[i] = (STACK *)createStack(100,sizeof(CTGinSCAF));

}
static void freeStackBuf(STACK **ctgStackBuffer,int scafBufSize)
{
	int i;
	for(i=0;i<scafBufSize;i++)
		freeStack(ctgStackBuffer[i]);
}

static void mapCtg2Scaf(int scafInBuf)
{
	int i,scafID;
	CTGinSCAF *actg;
	STACK *ctgsStack;
	unsigned int ctg,bal_ctg;
	
	for(i=0;i<scafInBuf;i++){
		scafID = scafCounter+i+1;
		ctgsStack = ctgStackBuffer[i];
		while((actg=stackPop(ctgsStack))!=NULL){
			ctg = actg->ctgID;
			bal_ctg = getTwinCtg(ctg);
			
			if(contig_array[ctg].from_vt!=0){
				contig_array[ctg].multi = 1;
				contig_array[bal_ctg].multi = 1;
				continue;
			}

			contig_array[ctg].from_vt = scafID;
			contig_array[ctg].to_vt = actg->start;
			contig_array[ctg].flag = 0;   //ctg and scaf on the same strand
			contig_array[bal_ctg].from_vt = scafID;
			contig_array[bal_ctg].to_vt = actg->start;
			contig_array[bal_ctg].flag = 1;
		}
	}

}

static void locateContigOnscaff(char *graphfile)
{

	FILE *fp;
	char line[1024];
	CTGinSCAF *actg;
	STACK *ctgStack,*aStack;
	int index=0,counter,overallLen;
	int starter,prev_start,gapN,scafLen;
	unsigned int ctg,prev_ctg=0;

	for(ctg=1;ctg<=num_ctg;ctg++){
		contig_array[ctg].from_vt = 0;
		contig_array[ctg].multi = 0;
	}

	ctgStack = (STACK *)createStack(1000,sizeof(CTGinSCAF));

	sprintf(line, "%s.scaf_gap", graphfile);
	fp = ckopen(line, "r");

	ctgStackBuffer = (STACK **)ckalloc(scafBufSize*sizeof(STACK *));
	initStackBuf(ctgStackBuffer,scafBufSize);


	Ncounter = scafCounter = scafInBuf = allGaps = 0;
	while(fgets(line,sizeof(line),fp)!=NULL){
		if(line[0]=='>'){
			if(index){
				aStack = ctgStackBuffer[scafInBuf++];
				reverseStack(aStack,ctgStack);
				if(scafInBuf==scafBufSize){
					mapCtg2Scaf(scafInBuf);
					scafCounter += scafInBuf;
					scafInBuf = 0;
				}
				//if(index%1000==0)
					//printf("Processed %d scaffolds\n",index);
			}
			//read next scaff	
			scafLen = prev_ctg = 0;
			emptyStack(ctgStack);
			sscanf(line+9,"%d %d %d",&index,&counter,&overallLen);
			fprintf(stderr,">%d\n",index);
			continue;
		}
		if(line[0]=='G'){  // gap appears
			continue;
		}
		if(line[0]>='0'&&line[0]<='9'){   // a contig line
			sscanf(line,"%d %d",&ctg,&starter);
			actg = (CTGinSCAF *)stackPush(ctgStack);	
			actg->ctgID = ctg;
			if(!prev_ctg){
				actg->start = scafLen;
				actg->end = actg->start + overlaplen + contig_array[ctg].length - 1;
			}else{
				gapN = starter - prev_start-(int)contig_array[prev_ctg].length;
				gapN = gapN < 1 ? 1:gapN;
				actg->start = scafLen + gapN;
				actg->end = actg->start + contig_array[ctg].length - 1;
			}
			fprintf(stderr,"%d\t%d\n",actg->start,actg->end);
			scafLen = actg->end+1;
			prev_ctg = ctg;
			prev_start = starter;
		}
	}
	if(index){
		aStack = ctgStackBuffer[scafInBuf++];
		reverseStack(aStack,ctgStack);
		mapCtg2Scaf(scafInBuf);
	}
	gapN = 0;
	for(ctg=1;ctg<=num_ctg;ctg++){
		if(contig_array[ctg].from_vt==0||contig_array[ctg].multi==1)
			continue;
		gapN++;
	}
	//printf("\nDone with %d scaffolds, %d contigs in Scaffolld\n",index,gapN);
	fclose(fp);
	freeStack(ctgStack);
	freeStackBuf(ctgStackBuffer,scafBufSize);
	free((void*)ctgStackBuffer);
}

static boolean contigElligible(unsigned int contigno)
{
	unsigned int ctg = index_array[contigno];
	if(contig_array[ctg].from_vt==0||contig_array[ctg].multi==1)
		return 0;	
	else
		return 1;	

}
static void output1read(FILE *fo,long long readno,unsigned int contigno,int pos)
{
	
	unsigned int ctg = index_array[contigno];
	int posOnScaf;
	char orien;
	pos = pos < 0 ? 0:pos;
	if(contig_array[ctg].flag==0){
		posOnScaf = contig_array[ctg].to_vt + pos - overlaplen;
		orien = '+';
	}else{
		posOnScaf = contig_array[ctg].to_vt + contig_array[ctg].length - pos;
		orien = '-';
	}
	/*
	if(readno==676)
		printf("Read %lld in region from %d, extend %d, pos %d, orien %c\n",
			readno,contig_array[ctg].to_vt,contig_array[ctg].length,posOnScaf,orien);
	*/
	fprintf(fo,"%lld\t%d\t%d\t%c\n",readno,contig_array[ctg].from_vt,posOnScaf,orien);
}

void locateReadOnScaf(char *graphfile)
{
	char name[1024],line[1024];
	FILE *fp,*fo;
	long long readno,counter=0,pre_readno=0;
	unsigned int contigno,pre_contigno;
	int pre_pos,pos;

	locateContigOnscaff(graphfile);

	sprintf(name,"%s.readOnContig",graphfile);
	fp = ckopen(name,"r");
	sprintf(name,"%s.readOnScaf",graphfile);
	fo = ckopen(name,"w");

	if(!orig2new){
		convertIndex();
		orig2new = 1;
	}
	fgets(line,1024,fp);
	while(fgets(line,1024,fp)!=NULL){
		sscanf(line,"%lld %d %d",&readno,&contigno,&pos);
		if((readno%2==0)&&(pre_readno==readno-1) // they are a pair of reads
			&&contigElligible(pre_contigno)&&contigElligible(contigno)){
			output1read(fo,pre_readno,pre_contigno,pre_pos);
			output1read(fo,readno,contigno,pos);
			counter++;
		}
		pre_readno = readno;
		pre_contigno = contigno;
		pre_pos = pos;
	}
	printf("%lld pairs on contig\n",counter);
	fclose(fp);
	fclose(fo);
}

#include "stdinc.h"
#include "newhash.h"
#include "extfunc.h"
#include "extvab.h"

#define CTGendLen 35  // shouldn't larger than max_read_len
#define UPlimit 5000
#define MaxRouteNum 10

static Kmer pubKmer = 0x1b4d65165b;

static void kmerSet_mark(KmerSet *set);
static void trace4Repeat(Kmer currW,int steps,int min,int max,int *num_route,
			KmerSet *kset,Kmer kmerDest,int overlap,Kmer WORDF,
			int *traceCounter,int maxRoute,kmer_t **soFarNode,short *multiOccu1,short *multiOccu2,
				int *routeLens,char **foundRoutes,char *soFarSeq,
			long long *soFarLinks,double *avgLinks);

static Kmer prevKmerLocal(Kmer next,char ch,int overlap)
{
	Kmer word = next;
	word >>= 2;
	word += ((Kmer)ch) << 2*(overlap-1);
	return word;
}
static Kmer nextKmerLocal(Kmer prev,char ch,Kmer WordFilter)
{
	Kmer word = prev;
	word <<= 2;
	word &= WordFilter;
	word += ch;
	return word;
}
static void singleKmer(int t, KmerSet *kset,int flag,Kmer *kmerBuffer,char *prevcBuffer,char *nextcBuffer)
{
	kmer_t *pos;
		
	put_kmerset(kset, kmerBuffer[t], prevcBuffer[t],nextcBuffer[t],&pos);
	if(pos->inEdge==flag)
		return;
	else if(pos->inEdge==0)
		pos->inEdge = flag;
	else if(pos->inEdge==1&&flag==2)
		pos->inEdge = 3;
	else if(pos->inEdge==2&&flag==1)
		pos->inEdge = 3;
	
}

static void putKmer2DBgraph(KmerSet *kset,int flag,int kmer_c,Kmer *kmerBuffer,char *prevcBuffer,char *nextcBuffer)
{
	int t;
	for(t=0;t<kmer_c;t++)
		singleKmer(t,kset,flag,kmerBuffer,prevcBuffer,nextcBuffer);

}

static void getSeqFromRead(READNEARBY read,char *src_seq)
{
	int len_seq=read.len;
	int j;
	char *tightSeq = (char *)darrayGet(readSeqInGap,read.seqStarter);
	for(j=0;j<len_seq;j++)
		src_seq[j] = getCharInTightString(tightSeq,j);
}

static void chopKmer4Ctg(Kmer *kmerCtg,int lenCtg,int overlap,char *src_seq,Kmer WORDF)
{
	int index,j;
	Kmer word = 0;
	for (index = 0;index<overlap;index++){
		word <<= 2;
		word += src_seq[index];
	}
	index = 0;
	kmerCtg[index++] = word;
	for(j = 1; j <= lenCtg - overlap; j ++)	{
		word = nextKmerLocal(word,src_seq[j-1+overlap],WORDF);
		kmerCtg[index++] = word;
	}
}

static void chopKmer4read(int len_seq,int overlap,char *src_seq,char *bal_seq,
		Kmer *kmerBuffer,char *prevcBuffer,char *nextcBuffer,int *kmer_c,Kmer WORDF)
{
	int j,bal_j;
	Kmer word,bal_word;
	int index;
	char InvalidCh=4;
	
	if(len_seq<overlap+1){
		*kmer_c = 0;
		return;
	}
	word = 0;
	for (index = 0;index<overlap;index++){
		word <<= 2;
		word += src_seq[index];
	}
	reverseComplementSeq(src_seq, len_seq,bal_seq);

		// complementary node
	bal_word = reverseComplement(word,overlap);
	bal_j = len_seq-0-overlap;  //  0;     
	index = 0;
	if(word<bal_word){
		kmerBuffer[index] = word;
		prevcBuffer[index] = InvalidCh;
		nextcBuffer[index++] = src_seq[0+overlap];
	}else{
		kmerBuffer[index] = bal_word;
		prevcBuffer[index] = bal_seq[bal_j-1];
		nextcBuffer[index++] = InvalidCh;
	}

	for(j = 1; j <= len_seq - overlap; j ++)	{
		word = nextKmerLocal(word,src_seq[j-1+overlap],WORDF);
		bal_j = len_seq-j-overlap; //  j; 
		bal_word = prevKmerLocal(bal_word,bal_seq[bal_j],overlap);
		
		if(word<bal_word){
			kmerBuffer[index] = word;
			prevcBuffer[index] = src_seq[j-1];
			if(j<len_seq - overlap)
				nextcBuffer[index++] = src_seq[j+overlap];
			else
				nextcBuffer[index++] = InvalidCh;
				//printf("%dth: %p with %p\n",kmer_c-1,word,hashBanBuffer[kmer_c-1]);
		}else{
				// complementary node
			kmerBuffer[index] = bal_word;
			if(bal_j>0)
				prevcBuffer[index] = bal_seq[bal_j-1];
			else
				prevcBuffer[index] = InvalidCh;
			nextcBuffer[index++] = bal_seq[bal_j+overlap];
				//printf("%dth: %p with %p\n",kmer_c-1,bal_word,hashBanBuffer[kmer_c-1]);
		}
	}
	*kmer_c = index;
}

static void headTightStr(char *tightStr,int length,int start,int headLen,int revS,char *src_seq)
{
	int i,index=0;

	if(!revS){
		for(i=start;i<start+headLen;i++)
			src_seq[index++] = getCharInTightString(tightStr,i);
	}
	else{
		for(i=length-1-start;i>=length-headLen-start;i--)
			src_seq[index++] = int_comp(getCharInTightString(tightStr,i));
	}
}

static int getSeqFromCtg(CTGinSCAF *ctg,boolean fromHead,unsigned int len,int originOverlap,char *src_seq)
{
	unsigned int ctgId = ctg->ctgID;
	unsigned int bal_ctg = getTwinCtg(ctgId);

	if(contig_array[ctgId].length<1)
		return 0;
	unsigned int length = contig_array[ctgId].length + originOverlap;
	
	len = len < length ? len:length;
	if(fromHead){
		if(contig_array[ctgId].seq)
			headTightStr(contig_array[ctgId].seq,length,0,len,0,src_seq);
		else
			headTightStr(contig_array[bal_ctg].seq,length,0,len,1,src_seq);
	}else{
		if(contig_array[ctgId].seq)
			headTightStr(contig_array[ctgId].seq,length,length-len,len,0,src_seq);
		else
			headTightStr(contig_array[bal_ctg].seq,length,length-len,len,1,src_seq);
	}
	return len;
}

			
static KmerSet *readsInGap2DBgraph(READNEARBY *rdArray, int num, CTGinSCAF *ctg1,CTGinSCAF *ctg2,int originOverlap,
				Kmer *kmerCtg1,Kmer *kmerCtg2,int overlap,Kmer WordFilter)
{
	int kmer_c;
	Kmer *kmerBuffer;
	char *nextcBuffer,*prevcBuffer;
	int i;
	int buffer_size=maxReadLen > CTGendLen ? maxReadLen:CTGendLen;
	KmerSet *kmerS=NULL;
	int lenCtg1;
	int lenCtg2;
	char *bal_seq;
	char *src_seq;

	src_seq = (char *)ckalloc(buffer_size*sizeof(char));
	bal_seq = (char *)ckalloc(buffer_size*sizeof(char));

	kmerBuffer = (Kmer *)ckalloc(buffer_size*sizeof(Kmer));
	prevcBuffer = (char *)ckalloc(buffer_size*sizeof(char));
	nextcBuffer = (char *)ckalloc(buffer_size*sizeof(char));

	kmerS = init_kmerset(1024,0.77f);

	for(i=0;i<num;i++){
		getSeqFromRead(rdArray[i],src_seq);
		chopKmer4read(rdArray[i].len,overlap,src_seq,bal_seq,
					kmerBuffer,prevcBuffer,nextcBuffer,&kmer_c,WordFilter);
		putKmer2DBgraph(kmerS,0,kmer_c,kmerBuffer,prevcBuffer,nextcBuffer);
	}

	lenCtg1 = getSeqFromCtg(ctg1,0,CTGendLen,originOverlap,src_seq);
	chopKmer4Ctg(kmerCtg1,lenCtg1,overlap,src_seq,WordFilter);
	chopKmer4read(lenCtg1,overlap,src_seq,bal_seq,
					kmerBuffer,prevcBuffer,nextcBuffer,&kmer_c,WordFilter);
	putKmer2DBgraph(kmerS,1,kmer_c,kmerBuffer,prevcBuffer,nextcBuffer);

	lenCtg2 = getSeqFromCtg(ctg2,1,CTGendLen,originOverlap,src_seq);
	chopKmer4Ctg(kmerCtg2,lenCtg2,overlap,src_seq,WordFilter);
	chopKmer4read(lenCtg2,overlap,src_seq,bal_seq,
					kmerBuffer,prevcBuffer,nextcBuffer,&kmer_c,WordFilter);
	putKmer2DBgraph(kmerS,2,kmer_c,kmerBuffer,prevcBuffer,nextcBuffer);
/*
	if(ctg1->ctgID==3733&&ctg2->ctgID==3067){
		for(i=0;i<lenCtg2;i++)
			printf("%c",int2base((int)src_seq[i]));
		printf("\n");
	}
*/
	//printf("sequence length chop from contigs on both sides: %d %d\n",lenCtg1,lenCtg2);
	//kmerSet_deLoop(kmerS,WordFilter);
	kmerSet_mark(kmerS);
	free((void *)src_seq);
	free((void *)bal_seq);
	free((void *)kmerBuffer);
	free((void *)nextcBuffer);
	free((void *)prevcBuffer);
	
	fflush(stdout);

	return kmerS;
}

static void printKmer(FILE *fo,Kmer kmer,int overlap)
{
	int i;
	char kmerSeq[32],ch;
	for(i=overlap-1;i>=0;i--){
		ch = kmer&3;
		kmer >>= 2;
		kmerSeq[i] = ch;
	}
	for(i=0;i<overlap;i++)
		fprintf(fo,"%c",int2base((int)kmerSeq[i]));
}

static void kmerSet_mark(KmerSet *set)
{
	int i,in_num,out_num,cvgSingle;
	kmer_t *rs;
	long long counter = 0,linear=0;
	Kmer word;

	set->iter_ptr = 0;
	while(set->iter_ptr < set->size){
		if(!is_kmer_entity_null(set->flags, set->iter_ptr)){
			in_num = out_num = 0;
			rs = set->array + set->iter_ptr;
			word = rs->seq;
			for(i=0;i<4;i++){
				cvgSingle = get_kmer_left_cov(*rs,i);
				if(cvgSingle>0){
					in_num++;
				}
				cvgSingle = get_kmer_right_cov(*rs,i);
				if(cvgSingle>0){
					out_num++;
				}
			}

			if(rs->single){
				counter++;
			}
			if(in_num==1&&out_num==1){
				rs->linear = 1;
				linear++;
			}
		}
		set->iter_ptr ++;
	}
	//printf("Allocated %ld node, %ld single nodes, %ld linear\n",(long)count_kmerset(set),counter,linear);
}

static kmer_t *searchNode(Kmer word,KmerSet *kset,int overlap)
{
	Kmer bal_word = reverseComplement(word,overlap);
	kmer_t *node;
	boolean found;
	if(word<bal_word)
		found = search_kmerset(kset,word,&node);	
	else
		found = search_kmerset(kset,bal_word,&node);	
	if(found)
		return node;
	else
		return NULL;
}

static int searchKmerOnCtg(Kmer currW,Kmer *kmerDest,int num)
{
	int i;
	for(i=0;i<num;i++){
		if(currW==kmerDest[i]){
			return i;
		}
	}
	return -1;
}

// pick on from n items randomly
static int nPick1(int *array,int n)
{
	int m,i;
	m = n-1;//(int)(drand48()*n);
	int value = array[m];
	for(i=m;i<n-1;i++)
		array[i] = array[i+1];
	return value;
}

static void traceAlongDBgraph(Kmer currW,int steps,int min,int max,int *num_route,
			KmerSet *kset,Kmer *kmerDest,int num,int overlap,Kmer WORDF,
			char **foundRoutes,int *routeEndOnCtg2,int *routeLens,char *soFarSeq,
			int *traceCounter,int maxRoute,kmer_t **soFarNode,boolean *multiOccu,
			long long *soFarLinks,double *avgLinks)
{
	(*traceCounter)++;
	if(*traceCounter>UPlimit){
		/*
		if(overlap==19&&kmerDest[0]==pubKmer)
			printf("UPlimit\n");
		*/
		return;
	}
	if(steps>max||*num_route>=maxRoute){
		/*
		if(overlap==19&&kmerDest[0]==pubKmer)
			printf("max steps/maxRoute\n");
		*/
		return;
	}
	Kmer word = reverseComplement(currW,overlap);
	boolean isSmaller = currW < word;
	int i;
	char ch;
	unsigned char links;
	if(isSmaller)
		word = currW;
	
	kmer_t *node;
	boolean found = search_kmerset(kset,word,&node);
	if(!found){
		printf("Trace: can't find kmer %llx (rc %llx, input %llx) at step %d\n",word,
				reverseComplement(word,overlap),currW,steps);
		return;
	}
	
	if(node->twin>1)
		return;
	if(soFarNode)
		soFarNode[steps] = node;

	if(steps>0)
		soFarSeq[steps-1] = currW&0x03;

	int index,end;
	int linkCounter = *soFarLinks;
	if(steps>=min&&node->inEdge>1&&(end=searchKmerOnCtg(currW,kmerDest,num))>=0){
		index = *num_route; 
		if(steps>0)
			avgLinks[index] = (double)linkCounter/steps;
		else
			avgLinks[index] = 0;
		//find node that appears more than once in the path
		multiOccu[index] = 0;
		for(i=0;i<steps+1;i++)
			soFarNode[i]->deleted = 0;
		for(i=0;i<steps+1;i++){
			if(soFarNode[i]->deleted){
				multiOccu[index] = 1;
				break;
			}
			soFarNode[i]->deleted = 1;
		}
		
		routeEndOnCtg2[index] = end;
		routeLens[index] = steps;
		char *array = foundRoutes[index]; 
		for(i=0;i<steps;i++)
			array[i] = soFarSeq[i];
		if(i<max)
			array[i] = 4;  //indicate the end of the sequence
		*num_route = ++index;
		return;
	}
	steps++;

	if(isSmaller){
		int array[] = {0,1,2,3};
		for(i=4;i>0;i--){
			ch = nPick1(array,i);
			links = get_kmer_right_cov(*node,ch);
			if(!links)
				continue;
			*soFarLinks = linkCounter + links;
			word = nextKmerLocal(currW,ch,WORDF);
			traceAlongDBgraph(word,steps,min,max,num_route,
						kset,kmerDest,num,overlap,WORDF,
					foundRoutes,routeEndOnCtg2,routeLens,soFarSeq,
						traceCounter,maxRoute,soFarNode,multiOccu,
					soFarLinks,avgLinks);
		}
	}else{
		int array[] = {0,1,2,3};
		for(i=4;i>0;i--){
                        ch = nPick1(array,i);
			links = get_kmer_left_cov(*node,ch);
			if(!links)
				continue;
			*soFarLinks = linkCounter + links;
			word = nextKmerLocal(currW,int_comp(ch),WORDF);
			traceAlongDBgraph(word,steps,min,max,num_route,
						kset,kmerDest,num,overlap,WORDF,
					foundRoutes,routeEndOnCtg2,routeLens,soFarSeq,
						traceCounter,maxRoute,soFarNode,multiOccu,
					soFarLinks,avgLinks);
		}
	}
}

static int searchFgap(KmerSet *kset,CTGinSCAF *ctg1,CTGinSCAF *ctg2,Kmer *kmerCtg1,
					Kmer *kmerCtg2,unsigned int origOverlap,int overlap,DARRAY *gapSeqArray,
			int len1,int len2,Kmer WordFilter,int *offset1,int *offset2,char *seqGap,int *cut1,int *cut2)
{

	int i;
	int ret = 0;
	kmer_t *node,**soFarNode;
	int num_route;
	int gapLen = ctg2->start - ctg1->end - origOverlap + overlap;
	int min = gapLen-GLDiff>0 ? gapLen-GLDiff:0; //0531
	int max = gapLen + GLDiff < 10 ? 10 : gapLen + GLDiff;
	char **foundRoutes;
	char *soFarSeq;
	int traceCounter;
	int *routeEndOnCtg2;
	int *routeLens;
	boolean *multiOccu;
	long long soFarLinks;
	double *avgLinks;

	//mask linear internal linear kmer on contig1 end
	routeEndOnCtg2 = (int *)ckalloc(MaxRouteNum*sizeof(int));
	routeLens = (int *)ckalloc(MaxRouteNum*sizeof(int));
	multiOccu = (boolean *)ckalloc(MaxRouteNum*sizeof(boolean));
	short *MULTI1 = (short *)ckalloc(MaxRouteNum*sizeof(short));
	short *MULTI2 = (short *)ckalloc(MaxRouteNum*sizeof(short));
	soFarSeq = (char *)ckalloc(max*sizeof(char));
	soFarNode = (kmer_t **)ckalloc((max+1)*sizeof(kmer_t *));
	foundRoutes = (char **)ckalloc(MaxRouteNum*sizeof(char *));;
	avgLinks = (double *)ckalloc(MaxRouteNum*sizeof(double));;
	for(i=0;i<MaxRouteNum;i++)
		foundRoutes[i] = (char *)ckalloc(max*sizeof(char));
	for(i=len1-1;i>=0;i--){
		
		num_route = traceCounter = soFarLinks = 0;
		int steps=0;
		traceAlongDBgraph(kmerCtg1[i],steps,min,max,&num_route,
				kset,kmerCtg2,len2,overlap,WordFilter,
					foundRoutes,routeEndOnCtg2,routeLens,soFarSeq,
					&traceCounter,MaxRouteNum,soFarNode,multiOccu,
				&soFarLinks,avgLinks);
		if(num_route>0){
			int m,minEnd=routeEndOnCtg2[0];
			for(m=0;m<num_route;m++){
				if(routeLens[m]<0)
					continue;
				if(routeEndOnCtg2[m]<minEnd)
					minEnd = routeEndOnCtg2[m];
			}
			 /* else if(minFreq>1){
				for(m=0;m<num_route;m++){
					if(routeEndOnCtg2[m]!=minEnd)
						continue;
					for(j=0;j<max;j++){
						if(foundRoutes[m][j]>3)
							break;
						printf("%c",int2base((int)foundRoutes[m][j]));
					}
					printf(": %4.2f\n",avgLinks[m]);
				}
			}   */
			
			num_route = traceCounter = soFarLinks = 0;
			steps=0;
			trace4Repeat(kmerCtg1[i],steps,min,max,&num_route,
				kset,kmerCtg2[minEnd],overlap,WordFilter,
					&traceCounter,MaxRouteNum,soFarNode,MULTI1,MULTI2,
				routeLens,foundRoutes,soFarSeq,&soFarLinks,avgLinks);
			int j,best=0;
			int maxLen=routeLens[0];
			double maxLink = avgLinks[0];
			char *pt;
			boolean repeat=0,sameLen=1;
			int leftMost=max,rightMost=max;
			if(num_route<1){
				fprintf(stderr,"After trace4Repeat: non route was found\n");
				continue;
			}
			if(num_route>1){
			// if multi paths are found, we check on the repeatative occurrences and links/length
				for(m=0;m<num_route;m++){
					if(routeLens[m]<0)
						continue;
					if(MULTI1[m]>=0&&MULTI2[m]>=0){
						repeat = 1;
						leftMost = leftMost>MULTI1[m] ? MULTI1[m]:leftMost;
						rightMost = rightMost>MULTI2[m] ? MULTI2[m]:rightMost;
					}
					if(routeLens[m]!=maxLen)
						sameLen = 0;
					if(routeLens[m]<maxLen)
						maxLen = routeLens[m];
					if(avgLinks[m]>maxLink){
						maxLink = avgLinks[m];
						best = m;
					}
				}
			}

			if(repeat){
				*offset1 = *offset2 = *cut1 = *cut2 = 0;
				int index=0;
				char ch;
				for(j=0;j<leftMost;j++){
					if(routeLens[0]<j+overlap+1)
						break;
					else
						ch=foundRoutes[0][j];
					for(m=1;m<num_route;m++){
						if(routeLens[m]<0)
							continue;
						if(ch!=foundRoutes[m][j])
							break;
					}
					if(m==num_route)
						seqGap[index++] = ch;
					else break;
				}
				
				*offset1 = index;
				index = 0;
				for(j=0;j<rightMost;j++){
					if(routeLens[0]-overlap-1<j)
						break;
					else
						ch=foundRoutes[0][routeLens[0]-overlap-1-j];
					for(m=1;m<num_route;m++){
						if(routeLens[m]<0)
							continue;
						if(ch!=foundRoutes[m][routeLens[m]-overlap-1-j])
							break;
					}
					if(m==num_route)
						index++;
					else break;
				}
				*offset2 = index;
				for(j=0;j<*offset2;j++)
					seqGap[*offset1+*offset2-1-j] = foundRoutes[0][routeLens[0]-overlap-1-j];
				if(*offset1>0||*offset2>0){
					*cut1 = len1-i-1;
					*cut2 = minEnd;
					//fprintf(stderr,"\n");
					for(m=0;m<num_route;m++){
						for(j=0;j<max;j++){
							if(foundRoutes[m][j]>3)
								break;
							//fprintf(stderr,"%c",int2base((int)foundRoutes[m][j]));
						}
						//fprintf(stderr,": %4.2f\n",avgLinks[m]);
					}
					/*
					fprintf(stderr,">Gap (%d + %d) (%d + %d)\n",*offset1,*offset2,*cut1,*cut2);
					for(index=0;index<*offset1+*offset2;index++)
						fprintf(stderr,"%c",int2base(seqGap[index]));
					fprintf(stderr,"\n");  */
				}
			
				ret = 3;
				break;
			}

			if(overlap+(len1-i-1)+minEnd-routeLens[best]>(int)origOverlap)
				continue;

			ctg1->gapSeqOffset = gapSeqArray->item_c;
			ctg1->gapSeqLen = routeLens[best];
			if(!darrayPut(gapSeqArray,ctg1->gapSeqOffset+maxLen/4))
				continue;
			pt = (char *)darrayPut(gapSeqArray,ctg1->gapSeqOffset);
			/*
			printKmer(stderr,kmerCtg1[i],overlap);
			fprintf(stderr,"-");
			*/
			for(j=0;j<max;j++){
				if(foundRoutes[best][j]>3)
					break;
				writeChar2tightString(foundRoutes[best][j],pt,j);
				//fprintf(stderr,"%c",int2base((int)foundRoutes[best][j]));
			}
			//fprintf(stderr,": GAPSEQ %d + %d, avglink %4.2f\n",len1-i-1,minEnd,avgLinks[best]);
			ctg1->cutTail = len1-i-1;
			ctg2->cutHead = overlap + minEnd;
			ctg2->scaftig_start = 0;

			ret = 1;
			break;
		/* }if(num_route>1){
			ret = 2;
			break;    */
		}else{  //mark node which leads to dead end
			node = searchNode(kmerCtg1[i],kset,overlap);
			if(node)
				node->twin = 2;  
		}
		
	}
	for(i=0;i<MaxRouteNum;i++)
		free((void *)foundRoutes[i]);
	free((void *)soFarSeq);
	free((void *)soFarNode);
	free((void *)multiOccu);
	free((void *)MULTI1);
	free((void *)MULTI2);
	free((void *)foundRoutes);
	free((void *)routeEndOnCtg2);
	free((void *)routeLens);

	return ret;
}

static void trace4Repeat(Kmer currW,int steps,int min,int max,int *num_route,
			KmerSet *kset,Kmer kmerDest,int overlap,Kmer WORDF,
			int *traceCounter,int maxRoute,kmer_t **soFarNode,short *multiOccu1,short *multiOccu2,
				int *routeLens,char **foundRoutes,char *soFarSeq,
			long long *soFarLinks,double *avgLinks)
{
	(*traceCounter)++;
	if(*traceCounter>UPlimit)
		return;
	if(steps>max||*num_route>=maxRoute)
		return;
	Kmer word = reverseComplement(currW,overlap);
	boolean isSmaller = currW < word;
	char ch;
	unsigned char links;
	int index,i;

	if(isSmaller)
		word = currW;
	
	kmer_t *node;
	boolean found = search_kmerset(kset,word,&node);
	if(!found){
		printf("Trace: can't find kmer %llx (rc %llx, input %llx) at step %d\n",word,
				reverseComplement(word,overlap),currW,steps);
		return;
	}
	if(soFarNode)
		soFarNode[steps] = node;
	if(soFarSeq&&steps>0)
		soFarSeq[steps-1] = currW&0x03;
	int linkCounter;
	if(soFarLinks)
		linkCounter = *soFarLinks;
	if(steps>=min&&currW==kmerDest){
		index = *num_route; 
		if(avgLinks&&steps>0)
			avgLinks[index] = (double)linkCounter/steps;
		else if(avgLinks)
			avgLinks[index] = 0;
		//find node that appears more than once in the path
		if(multiOccu1&&multiOccu2){
			for(i=0;i<steps+1;i++)
				soFarNode[i]->deleted = 0;
			int rightMost=0;
			boolean MULTI=0;
			for(i=0;i<steps+1;i++){
				if(soFarNode[i]->deleted){
					rightMost = rightMost<i-1 ? i-1:rightMost;
					MULTI = 1;
				}
				soFarNode[i]->deleted = 1;
			}
			if(!MULTI)
				multiOccu1[index] = multiOccu2[index] = -1;
			else{
				multiOccu2[index] = steps-2-rightMost<0 ? 0:steps-2-rightMost; //[0 steps-2]
				for(i=0;i<steps+1;i++)
					soFarNode[i]->deleted = 0;
				int leftMost=steps-2;
				for(i=steps;i>=0;i--){
					if(soFarNode[i]->deleted)
						leftMost = leftMost>i-1 ? i-1:leftMost;
					soFarNode[i]->deleted = 1;
				}
				multiOccu1[index] = leftMost<0 ? 0:leftMost;   //[0 steps-2]
			}
		}
		if(routeLens)
			routeLens[index] = steps;
		if(soFarSeq){
			char *array = foundRoutes[index]; 
			for(i=0;i<steps;i++)
				array[i] = soFarSeq[i];
			if(i<max)
				array[i] = 4;  //indicate the end of the sequence
		}
		*num_route = ++index;
	}

	steps++;

	if(isSmaller){
		int array[] = {0,1,2,3};
		for(i=4;i>0;i--){
                        ch = nPick1(array,i);
			links = get_kmer_right_cov(*node,ch);
			if(!links)
				continue;
			if(soFarLinks)
				*soFarLinks = linkCounter + links;
			word = nextKmerLocal(currW,ch,WORDF);
			trace4Repeat(word,steps,min,max,num_route,
						kset,kmerDest,overlap,WORDF,traceCounter,maxRoute,soFarNode,
					multiOccu1,multiOccu2,routeLens,foundRoutes,soFarSeq,
					soFarLinks,avgLinks);
		}
	}else{
		int array[] = {0,1,2,3};
		for(i=4;i>0;i--){
                        ch = nPick1(array,i);
			links = get_kmer_left_cov(*node,ch);
			if(!links)
				continue;
			if(soFarLinks)
				*soFarLinks = linkCounter + links;
			word = nextKmerLocal(currW,int_comp(ch),WORDF);
			trace4Repeat(word,steps,min,max,num_route,
						kset,kmerDest,overlap,WORDF,traceCounter,maxRoute,soFarNode,
					multiOccu1,multiOccu2,routeLens,foundRoutes,soFarSeq,
					soFarLinks,avgLinks);
		}
	}
}

//found repeat node on contig ends
static void maskRepeatNode(KmerSet *kset,Kmer *kmerCtg1,
					Kmer *kmerCtg2,int overlap,
					int len1,int len2,int max,Kmer WordFilter)
{
	int i;
	int num_route,steps;
	int min = 1,maxRoute=1;
	int traceCounter;
	Kmer word,bal_word;
	kmer_t *node;
	boolean found;
	int counter=0;
	for(i=0;i<len1;i++){
		word = kmerCtg1[i];	
		bal_word = reverseComplement(word,overlap);
		if(word>bal_word)
			word=bal_word;
		found = search_kmerset(kset,word,&node);	
		if(!found||node->linear){
			//printf("Found no node for kmer %llx\n",word);
			continue;
		}
		num_route = traceCounter = 0;
		steps=0;
		trace4Repeat(word,steps,min,max,&num_route,
				kset,word,overlap,WordFilter,
					&traceCounter,maxRoute,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
		if(num_route<1)
			continue;
		counter++;
		node->checked = 1;
	}
	for(i=0;i<len2;i++){
		word = kmerCtg2[i];	
		bal_word = reverseComplement(word,overlap);
		if(word>bal_word)
			word=bal_word;
		found = search_kmerset(kset,word,&node);	
		if(!found||node->linear){
			//printf("Found no node for kmer %llx\n",word);
			continue;
		}
		num_route = traceCounter = 0;
		steps=0;
		trace4Repeat(word,steps,min,max,&num_route,
				kset,word,overlap,WordFilter,
					&traceCounter,maxRoute,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
		if(num_route<1)
			continue;
		counter++;
		node->checked = 1;
	}
	//printf("MR: %d(%d)\n",counter,len1+len2);
}

/*
static boolean chopReadFillGap(int len_seq,int overlap,char *src_seq, char *bal_seq,
					KmerSet *kset,Kmer WORDF,int *start,int *end,boolean *bal, 
					Kmer *KmerCtg1,int len1,Kmer *KmerCtg2,int len2,int *index1,int *index2)
{
	int index,j=0,bal_j;
	Kmer word,bal_word;
	int flag=0,bal_flag=0;
	int ctg1start,bal_ctg1start,ctg2end,bal_ctg2end;
	int seqStart,bal_start,seqEnd,bal_end;
	kmer_t *node;
	boolean found;
	
	if(len_seq<overlap+1){
		return 0;
	}
	word = 0;
	for (index = 0;index<overlap;index++){
		word <<= 2;
		word += src_seq[index];
	}
	reverseComplementSeq(src_seq, len_seq,bal_seq);

		// complementary node
	bal_word = reverseComplement(word,overlap);
	bal_j = len_seq-0-overlap;  //  0;     
	flag = bal_flag = 0;
	if(word<bal_word){
		found = search_kmerset(kset,word,&node);	
	}else{
		found = search_kmerset(kset,bal_word,&node);	
	}
	if(found&&!node->linear&&!node->checked){
		if(!flag&&node->inEdge==1){
			ctg1start = searchKmerOnCtg(word,KmerCtg1,len1);
			if(ctg1start>0){
				flag = 1;	
				seqStart = j + overlap-1;
			}
		}	
		if(!bal_flag&&node->inEdge==2){
			bal_ctg2end = searchKmerOnCtg(bal_word,KmerCtg2,len2);
			if(bal_ctg2end>0){
				bal_flag = 2;	
				bal_end = bal_j+overlap-1;
			}
		}	
	}

	for(j = 1; j <= len_seq - overlap; j ++)	{
		word = nextKmerLocal(word,src_seq[j-1+overlap],WORDF);
		bal_j = len_seq-j-overlap; //  j; 
		bal_word = prevKmerLocal(bal_word,bal_seq[bal_j],overlap);
		
		if(word<bal_word){
			found = search_kmerset(kset,word,&node);	
		}else{
			found = search_kmerset(kset,bal_word,&node);	
		}
		if(found&&!node->linear&&!node->checked){
			if(!flag&&node->inEdge==1){
				ctg1start = searchKmerOnCtg(word,KmerCtg1,len1);
				if(ctg1start>0){
					flag = 1;	
					seqStart = j + overlap-1;
				}
			}else if(flag==1&&node->inEdge==1){
				index = searchKmerOnCtg(word,KmerCtg1,len1);
				if(index>ctg1start){    // choose hit closer to gap
					ctg1start = index;
					seqStart = j + overlap-1;
				}
			}else if(flag==1&&node->inEdge==2){
				ctg2end = searchKmerOnCtg(word,KmerCtg2,len2);
				if(ctg2end>0){
					flag = 3;	
					seqEnd = j+overlap-1;
					break;
				}
			}	

			if(!bal_flag&&node->inEdge==2){
				bal_ctg2end = searchKmerOnCtg(bal_word,KmerCtg2,len2);
				if(bal_ctg2end>0){
					bal_flag = 2;	
					bal_end = bal_j+overlap-1;
				}
			}else if(bal_flag==2&&node->inEdge==2){
				index = searchKmerOnCtg(bal_word,KmerCtg2,len2);
				if(index<bal_ctg2end){   // choose hit closer to gap
					index = bal_ctg2end;
					bal_end = bal_j+overlap-1;
				}
			}else if(bal_flag==2&&node->inEdge==1){
				bal_ctg1start = searchKmerOnCtg(bal_word,KmerCtg1,len1);
				if(bal_ctg1start>0){
					bal_flag = 3;	
					bal_start = bal_j+overlap-1;
					break;
				}
			}	
		}
	}
	if(flag==3){
		*start = seqStart;
		*end = seqEnd;
		*bal = 0;
		*index1 = ctg1start;
		*index2 = ctg2end;
		return 1;	
	}else if(bal_flag==3){
		*start = bal_start;
		*end = bal_end;
		*bal = 1;
		*index1 = bal_ctg1start;
		*index2 = bal_ctg2end;
		return 1;	
	}
	return 0;
}

static boolean readsCrossGap(READNEARBY *rdArray, int num, int originOverlap,DARRAY *gapSeqArray,
				Kmer *kmerCtg1,Kmer *kmerCtg2,int overlap,int len1,int len2,
					CTGinSCAF *ctg1,CTGinSCAF *ctg2,KmerSet *kmerS,Kmer WordFilter,int min,int max)
{
	int i,j,start,end,startOnCtg1,endOnCtg2;
	char *bal_seq;
	char *src_seq;
	char *pt;
	boolean bal,ret=0,FILL;
	
	src_seq = (char *)ckalloc(maxReadLen*sizeof(char));
	bal_seq = (char *)ckalloc(maxReadLen*sizeof(char));

	for(i=0;i<num;i++){
		getSeqFromRead(rdArray[i],src_seq);
		FILL = chopReadFillGap(rdArray[i].len,overlap,src_seq,bal_seq,
					kmerS,WordFilter,&start,&end,&bal,
					kmerCtg1,len1,kmerCtg2,len2,&startOnCtg1,&endOnCtg2);

		if(!FILL||(end-start)<min||(end-start)>max)
			continue;
		fprintf(stderr,"Read across\n");
		//printf("Filled: K %d, ctg1 %d ctg2 %d,start %d end %d\n",overlap,startOnCtg1,endOnCtg2,start,end);
		if(overlap+(len1-startOnCtg1-1)+endOnCtg2-(end-start)>(int)originOverlap)
			continue;  // contig1 and contig2 could not overlap more than origOverlap bases

		ctg1->gapSeqOffset = gapSeqArray->item_c;
		ctg1->gapSeqLen = end-start;
		if(!darrayPut(gapSeqArray,ctg1->gapSeqOffset+(end-start)/4))
			continue;
		pt = (char *)darrayPut(gapSeqArray,ctg1->gapSeqOffset);
		for(j=start+1;j<=end;j++){
			if(bal)
				writeChar2tightString(bal_seq[j],pt,j-start-1);
			else
				writeChar2tightString(src_seq[j],pt,j-start-1);
			
		}
		ctg1->cutTail = len1-startOnCtg1-1;
		ctg2->cutHead = overlap + endOnCtg2;
		ctg2->scaftig_start = 0;
		
		ret = 1;
		break;
	}

	free((void*)src_seq);
	free((void*)bal_seq);
	return ret;
}
*/
static void kmerSet_markTandem(KmerSet *set,Kmer WordFilter,int overlap);
static boolean readsCrossGap(READNEARBY *rdArray, int num, int originOverlap,DARRAY *gapSeqArray,
				Kmer *kmerCtg1,Kmer *kmerCtg2,int overlap,
					CTGinSCAF *ctg1,CTGinSCAF *ctg2,KmerSet *kmerS,Kmer WordFilter,int min,int max,
			int offset1,int offset2,char *seqGap,char *seqCtg1,char *seqCtg2,int cut1,int cut2);

int localGraph(READNEARBY *rdArray,int num,CTGinSCAF *ctg1,CTGinSCAF *ctg2,
						int origOverlap,Kmer *kmerCtg1,Kmer *kmerCtg2,
			int overlap,DARRAY *gapSeqArray,char *seqCtg1,char *seqCtg2,char *seqGap)
{
	/**************** put kmer in DBgraph ****************/
	KmerSet *kmerSet;
	Kmer WordFilter = (((Kmer) 1) << (2*overlap)) - 1;
/*
	if(ctg1->ctgID==56410&&ctg2->ctgID==61741)
		printf("Extract %d reads for gap [%d %d]\n",num,ctg1->ctgID,ctg2->ctgID);
*/
	kmerSet = readsInGap2DBgraph(rdArray,num,ctg1,ctg2,origOverlap,
					kmerCtg1,kmerCtg2,overlap,WordFilter);
	time_t tt;
	time(&tt);
//	srand48((int)tt);
/*
	int i,j;
	for(i=0;i<2;i++){
		int array[] = {0,1,2,3};
		for(j=4;j>0;j--)
			fprintf(stderr,"%d ", nPick1(array,j));
	}
	fprintf(stderr,"\n");
*/
	/***************** search path to connect contig ends ********/
	int gapLen = ctg2->start - ctg1->end - origOverlap + overlap;
	int min = gapLen-GLDiff>0 ? gapLen-GLDiff:0; 
	int max = gapLen + GLDiff < 10 ? 10 : gapLen + GLDiff;
	//count kmer number for contig1 and contig2 ends
	int len1,len2;
	len1 = CTGendLen<contig_array[ctg1->ctgID].length+origOverlap ? 
				CTGendLen:contig_array[ctg1->ctgID].length+origOverlap;
	len2 = CTGendLen<contig_array[ctg2->ctgID].length+origOverlap ? 
				CTGendLen:contig_array[ctg2->ctgID].length+origOverlap;
	len1 -= overlap-1;
	len2 -= overlap-1;
	
	//int pathNum = 2;
	int offset1=0,offset2=0,cut1=0,cut2=0;
	int pathNum = searchFgap(kmerSet,ctg1,ctg2,kmerCtg1,kmerCtg2,
				origOverlap,overlap,gapSeqArray,
				len1,len2,WordFilter,&offset1,&offset2,seqGap,&cut1,&cut2);
	
	//printf("SF: %d K %d\n",pathNum,overlap);
	if(pathNum==0){
		free_kmerset(kmerSet);
		return 0;
	}else if(pathNum==1){
		free_kmerset(kmerSet);
		return 1;
	}/*
	else{
		printf("ret %d\n",pathNum);
		free_kmerset(kmerSet);
		return 0;
	}  */

	/******************* cross the gap by single reads *********/
	//kmerSet_markTandem(kmerSet,WordFilter,overlap);
	maskRepeatNode(kmerSet,kmerCtg1,kmerCtg2,overlap,
					len1,len2,max,WordFilter);
	boolean found = readsCrossGap(rdArray,num,origOverlap,gapSeqArray,
				kmerCtg1,kmerCtg2,overlap,ctg1,ctg2,kmerSet,WordFilter,min,max,
			offset1,offset2,seqGap,seqCtg1,seqCtg2,cut1,cut2);
	if(found){
		//fprintf(stderr,"read across\n");
		free_kmerset(kmerSet);
		return found;
	}
	else{
		free_kmerset(kmerSet);
		return 0;
	}
	
}

static void kmerSet_markTandem(KmerSet *set,Kmer WordFilter,int overlap)
{
	kmer_t *rs;
	long long counter = 0;
	int num_route,steps;
	int min=1,max=overlap,maxRoute=1;
	int traceCounter;

	set->iter_ptr = 0;
	while(set->iter_ptr < set->size){
		if(!is_kmer_entity_null(set->flags, set->iter_ptr)){
			rs = set->array + set->iter_ptr;
			if(rs->inEdge>0){
				set->iter_ptr ++;
				continue;
			}
			num_route = traceCounter = 0;
			steps=0;
			trace4Repeat(rs->seq,steps,min,max,&num_route,
				set,rs->seq,overlap,WordFilter,
					&traceCounter,maxRoute,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
			if(num_route<1){
				set->iter_ptr ++;
				continue;
			}
			/*
			printKmer(stderr,rs->seq,overlap);
			fprintf(stderr, "\n");
			*/
			rs->checked = 1;
			counter++;
		}
		set->iter_ptr ++;
	}
}
/******************* the following is for read-crossing gaps *************************/

#define MAXREADLENGTH 100

static const int INDEL = 0;
static const int SIM[4][4] = {
	{1, 0, 0, 0},
	{0, 1, 0, 0},
	{0, 0, 1, 0},
	{0, 0, 0, 1}
};
static char fastSequence[MAXREADLENGTH];
static char slowSequence[MAXREADLENGTH];

static int Fmatrix[MAXREADLENGTH + 1][MAXREADLENGTH + 1];
static int slowToFastMapping[MAXREADLENGTH + 1];
static int fastToSlowMapping[MAXREADLENGTH + 1];

static int max(int A, int B, int C)
{
	A = A>=B ? A:B;
	return (A>=C ? A:C);

}

static int compareSequences(char * sequence1, char * sequence2, int length1, int length2)
{
	if(length1<1||length2<1||length1>MAXREADLENGTH||length2>MAXREADLENGTH)
		return 0;
	int i, j;
	int Choice1, Choice2, Choice3;
	int maxScore;

	for (i = 0; i <= length1; i++)
		Fmatrix[i][0] = 0;
	for (j = 0; j <= length2; j++)
		Fmatrix[0][j] = 0;

	for (i = 1; i <= length1; i++) {
		for (j = 1; j <= length2; j++) {
			Choice1 =
			    Fmatrix[i - 1][j - 1] +
			    SIM[(int) sequence1[i-1]]
			    [(int) sequence2[j-1]];
			Choice2 = Fmatrix[i - 1][j] + INDEL;
			Choice3 = Fmatrix[i][j - 1] + INDEL;
			Fmatrix[i][j] = max(Choice1, Choice2, Choice3);
		}
	}

	maxScore = Fmatrix[length1][length2];
	return maxScore;
}

static void mapSlowOntoFast(int slowSeqLength,int fastSeqLength)
{
	int slowIndex = slowSeqLength;
	int fastIndex = fastSeqLength;
	int fastn, slown;

	if (slowIndex == 0) {
		slowToFastMapping[0] = fastIndex;

		while (fastIndex >= 0)
			fastToSlowMapping[fastIndex--] = 0;

		return;
	}

	if (fastIndex == 0) {
		while (slowIndex >= 0)
			slowToFastMapping[slowIndex--] = 0;

		fastToSlowMapping[0] = slowIndex;

		return;
	}

	while (slowIndex > 0 && fastIndex > 0) {
		fastn = (int) fastSequence[fastIndex-1]; //getCharInTightString(fastSequence,fastIndex-1);
		slown = (int) slowSequence[slowIndex-1]; //getCharInTightString(slowSequence,slowIndex-1);

		if (Fmatrix[fastIndex][slowIndex] ==
		    Fmatrix[fastIndex - 1][slowIndex - 1] +
		    SIM[fastn][slown]) {
			fastToSlowMapping[--fastIndex] = --slowIndex;
			slowToFastMapping[slowIndex] = fastIndex;
		} else if (Fmatrix[fastIndex][slowIndex] ==
			   Fmatrix[fastIndex - 1][slowIndex] + INDEL)
			fastToSlowMapping[--fastIndex] = slowIndex - 1;

		else if (Fmatrix[fastIndex][slowIndex] ==
			 Fmatrix[fastIndex][slowIndex - 1] + INDEL)
			slowToFastMapping[--slowIndex] = fastIndex - 1;

		else {
			printf("compareSequence: Error trace\n");
			fflush(stdout);
			abort();
		}
	}

	while (slowIndex > 0)
		slowToFastMapping[--slowIndex] = -1;
	while (fastIndex > 0)
		fastToSlowMapping[--fastIndex] = -1;

	slowToFastMapping[slowSeqLength] =
	    fastSeqLength;
	fastToSlowMapping[fastSeqLength] =
	    slowSeqLength;
}

static boolean chopReadFillGap(int len_seq,int overlap,char *src_seq, char *bal_seq,
					KmerSet *kset,Kmer WORDF,int *start,int *end,boolean *bal, 
					Kmer *KmerCtg1,int len1,Kmer *KmerCtg2,int len2,int *index1,int *index2)
{
	int index,j=0,bal_j;
	Kmer word,bal_word;
	int flag=0,bal_flag=0;
	int ctg1start,bal_ctg1start,ctg2end,bal_ctg2end;
	int seqStart,bal_start,seqEnd,bal_end;
	kmer_t *node;
	boolean found;
	
	if(len_seq<overlap+1){
		return 0;
	}
	word = 0;
	for (index = 0;index<overlap;index++){
		word <<= 2;
		word += src_seq[index];
	}
	reverseComplementSeq(src_seq, len_seq,bal_seq);

		// complementary node
	bal_word = reverseComplement(word,overlap);
	bal_j = len_seq-0-overlap;  //  0;     
	flag = bal_flag = 0;
	if(word<bal_word){
		found = search_kmerset(kset,word,&node);	
	}else{
		found = search_kmerset(kset,bal_word,&node);	
	}
	if(found&&!node->linear&&!node->checked){
		if(!flag&&node->inEdge==1){
			ctg1start = searchKmerOnCtg(word,KmerCtg1,len1);
			if(ctg1start>=0){
				flag = 1;	
				seqStart = j + overlap-1;
			}
		}	
		if(!bal_flag&&node->inEdge==2){
			bal_ctg2end = searchKmerOnCtg(bal_word,KmerCtg2,len2);
			if(bal_ctg2end>=0){
				bal_flag = 2;	
				bal_end = bal_j+overlap-1;
			}
		}	
	}

	for(j = 1; j <= len_seq - overlap; j ++)	{
		word = nextKmerLocal(word,src_seq[j-1+overlap],WORDF);
		bal_j = len_seq-j-overlap; //  j; 
		bal_word = prevKmerLocal(bal_word,bal_seq[bal_j],overlap);
		
		if(word<bal_word){
			found = search_kmerset(kset,word,&node);	
		}else{
			found = search_kmerset(kset,bal_word,&node);	
		}
		if(found&&!node->linear&&!node->checked){
			if(!flag&&node->inEdge==1){
				ctg1start = searchKmerOnCtg(word,KmerCtg1,len1);
				if(ctg1start>=0){
					flag = 1;	
					seqStart = j + overlap-1;
				}
			}else if(flag==1&&node->inEdge==1){
				index = searchKmerOnCtg(word,KmerCtg1,len1);
				if(index>=0&&index>ctg1start){    // choose hit closer to gap
					ctg1start = index;
					seqStart = j + overlap-1;
				}
			}else if(flag==1&&node->inEdge==2){
				ctg2end = searchKmerOnCtg(word,KmerCtg2,len2);
				if(ctg2end>=0){
					flag = 3;	
					seqEnd = j+overlap-1;
					break;
				}
			}	

			if(!bal_flag&&node->inEdge==2){
				bal_ctg2end = searchKmerOnCtg(bal_word,KmerCtg2,len2);
				if(bal_ctg2end>=0){
					bal_flag = 2;	
					bal_end = bal_j+overlap-1;
				}
			}else if(bal_flag==2&&node->inEdge==2){
				index = searchKmerOnCtg(bal_word,KmerCtg2,len2);
				if(index>=0&&index<bal_ctg2end){   // choose hit closer to gap
					bal_ctg2end = index;
					bal_end = bal_j+overlap-1;
				}
			}else if(bal_flag==2&&node->inEdge==1){
				bal_ctg1start = searchKmerOnCtg(bal_word,KmerCtg1,len1);
				if(bal_ctg1start>=0){
					bal_flag = 3;	
					bal_start = bal_j+overlap-1;
					break;
				}
			}	
		}
	}
	if(flag==3){
		*start = seqStart;
		*end = seqEnd;
		*bal = 0;
		*index1 = ctg1start;
		*index2 = ctg2end;
		return 1;	
	}else if(bal_flag==3){
		*start = bal_start;
		*end = bal_end;
		*bal = 1;
		*index1 = bal_ctg1start;
		*index2 = bal_ctg2end;
		return 1;	
	}
	return 0;
}


static int cutSeqFromTightStr(char *tightStr,int length,int start,int end,int revS,char *src_seq)
{
	int i,index=0;
	end = end < length ? end:length-1;
	start = start>=0 ? start:0;

	if(!revS){
		for(i=start;i<=end;i++)
			src_seq[index++] = getCharInTightString(tightStr,i);
	}
	else{
		for(i=length-1-start;i>=length-end-1;i--)
			src_seq[index++] = int_comp(getCharInTightString(tightStr,i));
	}
	return end-start+1;
}

static int cutSeqFromCtg(unsigned int ctgID,int start,int end, char *sequence,int originOverlap)
{
	
	unsigned int bal_ctg = getTwinCtg(ctgID);
	if(contig_array[ctgID].length<1)
		return 0;
	int length = contig_array[ctgID].length+originOverlap;
	if(contig_array[ctgID].seq)
		return cutSeqFromTightStr(contig_array[ctgID].seq,length,start,end,0,sequence);
	else
		return cutSeqFromTightStr(contig_array[bal_ctg].seq,length,start,end,1,sequence);

}

static int cutSeqFromRead(char *src_seq,int length,int start,int end,char *sequence)
{
	if(end>=length)
		printf("******: end %d length %d\n",end,length);
	end = end<length ? end:length-1;
	start = start>=0 ? start:0;
	int i;
	for(i=start;i<=end;i++)
		sequence[i-start] = src_seq[i];
	return end-start+1;
}

static void printSeq(FILE *fo,char *seq,int len)
{
	int i;
	for(i=0;i<len;i++)
		fprintf(fo,"%c",int2base((int)seq[i]));

	fprintf(fo,"\n");
}

static boolean readsCrossGap(READNEARBY *rdArray, int num, int originOverlap,DARRAY *gapSeqArray,
				Kmer *kmerCtg1,Kmer *kmerCtg2,int overlap,
					CTGinSCAF *ctg1,CTGinSCAF *ctg2,KmerSet *kmerS,Kmer WordFilter,int min,int max,
			int offset1,int offset2,char *seqGap,char *seqCtg1,char *seqCtg2,int cut1,int cut2)
{
	int i,j,start,end,startOnCtg1,endOnCtg2;
	char *bal_seq;
	char *src_seq;
	char *pt;
	boolean bal,ret=0,FILL;
	double maxScore=0.0;
	int maxIndex;
	int lenCtg1,lenCtg2;
	//build sequences on left and right of the uncertain region 
	int buffer_size=maxReadLen > 100 ? maxReadLen:100;
	int length = contig_array[ctg1->ctgID].length+originOverlap;
	if(buffer_size>offset1){
		lenCtg1 = cutSeqFromCtg(ctg1->ctgID,length-cut1-(buffer_size-offset1),length-1-cut1,seqCtg1,originOverlap);
		for(i=0;i<offset1;i++)
			seqCtg1[lenCtg1+i] = seqGap[i];
		lenCtg1 += offset1;
	}else{
		for(i=offset1-buffer_size;i<offset1;i++)
			seqCtg1[i+buffer_size-offset1] = seqGap[i];
		lenCtg1 = buffer_size;
	}
	length = contig_array[ctg2->ctgID].length+originOverlap;
	if(buffer_size>offset2){
		lenCtg2 = cutSeqFromCtg(ctg2->ctgID,cut2,buffer_size-offset2-1+cut2,&(seqCtg2[offset2]),originOverlap);
		for(i=0;i<offset2;i++)
			seqCtg2[i] = seqGap[i+offset1];
		lenCtg2 += offset2;
	}else{
		for(i=0;i<buffer_size;i++)
			seqCtg2[i] = seqGap[i+offset1];
		lenCtg2 = buffer_size;
	}
		/*
	if(offset1>0||offset2>0){
		for(i=0;i<lenCtg1;i++)
			fprintf(stderr,"%c",int2base(seqCtg1[i]));
		fprintf(stderr,": CTG1\n");
		for(i=0;i<lenCtg2;i++)
			fprintf(stderr,"%c",int2base(seqCtg2[i]));
		fprintf(stderr,": CTG2\n");
	}
		*/
	//chop kmer from both ends of the uncertain region
	int len1,len2;
	len1 = CTGendLen<lenCtg1 ? CTGendLen:lenCtg1;
	len2 = CTGendLen<lenCtg2 ? CTGendLen:lenCtg2;
	chopKmer4Ctg(kmerCtg1,len1,overlap,&(seqCtg1[lenCtg1-len1]),WordFilter);
	chopKmer4Ctg(kmerCtg2,len2,overlap,seqCtg2,WordFilter);
	len1 -= overlap-1;
	len2 -= overlap-1;
	
	src_seq = (char *)ckalloc(maxReadLen*sizeof(char));
	bal_seq = (char *)ckalloc(maxReadLen*sizeof(char));
	
	int *START = (int *)ckalloc(num*sizeof(int));
	int *END = (int *)ckalloc(num*sizeof(int));
	int *INDEX1 = (int *)ckalloc(num*sizeof(int));
	int *INDEX2 = (int *)ckalloc(num*sizeof(int));
	double *SCORE = (double *)ckalloc(num*sizeof(double));
	boolean *BAL = (boolean *)ckalloc(num*sizeof(boolean));
	memset(SCORE,0,num*sizeof(double));
	
	for(i=0;i<num;i++){
		getSeqFromRead(rdArray[i],src_seq);
		FILL = chopReadFillGap(rdArray[i].len,overlap,src_seq,bal_seq,
					kmerS,WordFilter,&start,&end,&bal,
					kmerCtg1,len1,kmerCtg2,len2,&startOnCtg1,&endOnCtg2);

		if(!FILL||(end-start)<min||(end-start)>max)
			continue;
		if(overlap+(len1-startOnCtg1-1)+endOnCtg2-(end-start)>(int)originOverlap)
			continue;  // contig1 and contig2 could not overlap more than origOverlap bases
		START[i] = start;
		END[i] = end;
		INDEX1[i] = startOnCtg1;
		INDEX2[i] = endOnCtg2;
		BAL[i] = bal;

		int matchLen = 2*overlap<(end-start+overlap) ? 2*overlap:(end-start+overlap);
		int match;
		int alignLen = matchLen;
		//compare the left of hit kmer on ctg1
		//int ctgLeft = (contig_array[ctg1->ctgID].length+originOverlap)-(len1+overlap-1)+startOnCtg1;
		int ctgLeft = (lenCtg1)-(len1+overlap-1)+startOnCtg1;
		int readLeft = start-overlap+1;
		int cmpLen = ctgLeft<readLeft ? ctgLeft:readLeft;
		
		cmpLen = cmpLen<=MAXREADLENGTH ? cmpLen:MAXREADLENGTH;
		//cutSeqFromCtg(ctg1->ctgID,ctgLeft-cmpLen,ctgLeft-1,fastSequence,originOverlap);
		cutSeqFromRead(seqCtg1,lenCtg1,ctgLeft-cmpLen,ctgLeft-1,fastSequence);
		if(!bal)
			cutSeqFromRead(src_seq,rdArray[i].len,readLeft-cmpLen,readLeft-1,slowSequence);
		else
			cutSeqFromRead(bal_seq,rdArray[i].len,readLeft-cmpLen,readLeft-1,slowSequence);
		match = compareSequences(fastSequence,slowSequence, cmpLen, cmpLen);
		
		alignLen += cmpLen;
		matchLen += match;
		
		//compare the right of hit kmer on ctg1
		int ctgRight = len1-startOnCtg1-1;
		
		cmpLen = ctgRight<(rdArray[i].len-start-1) ? ctgRight:(rdArray[i].len-start-1);
		cmpLen = cmpLen<=MAXREADLENGTH ? cmpLen:MAXREADLENGTH;
		//cutSeqFromCtg(ctg1->ctgID,ctgLeft+overlap,ctgLeft+overlap+cmpLen-1,fastSequence,originOverlap);
		cutSeqFromRead(seqCtg1,lenCtg1,ctgLeft+overlap,ctgLeft+overlap+cmpLen-1,fastSequence);
		if(!bal)
			cutSeqFromRead(src_seq,rdArray[i].len,start+1,start+cmpLen,slowSequence);
		else
			cutSeqFromRead(bal_seq,rdArray[i].len,start+1,start+cmpLen,slowSequence);
		match = compareSequences(fastSequence,slowSequence, cmpLen, cmpLen);
		//fprintf(stderr,"%d -- %d\n",match,cmpLen);
		
		alignLen += cmpLen;
		matchLen += match;
		
		//compare the left of hit kmer on ctg2
		ctgLeft = endOnCtg2;
		readLeft = end-overlap+1;
		cmpLen = ctgLeft<readLeft ? ctgLeft:readLeft;
		cmpLen = ctgLeft<=MAXREADLENGTH ? ctgLeft:MAXREADLENGTH;
		//cutSeqFromCtg(ctg2->ctgID,endOnCtg2-cmpLen,endOnCtg2-1,fastSequence,originOverlap);
		cutSeqFromRead(seqCtg2,lenCtg2,endOnCtg2-cmpLen,endOnCtg2-1,fastSequence);
		if(!bal)
			cutSeqFromRead(src_seq,rdArray[i].len,readLeft-cmpLen,readLeft-1,slowSequence);
		else
			cutSeqFromRead(bal_seq,rdArray[i].len,readLeft-cmpLen,readLeft-1,slowSequence);
		match = compareSequences(fastSequence,slowSequence, cmpLen, cmpLen);
		alignLen += cmpLen;
		matchLen += match;
		
		//compare the right of hit kmer on ctg2
		//ctgRight = contig_array[ctg2->ctgID].length+originOverlap-endOnCtg2-overlap;
		ctgRight = lenCtg2-endOnCtg2-overlap;
		cmpLen = ctgRight<(rdArray[i].len-end-1) ? ctgRight:(rdArray[i].len-end-1);
		cmpLen = cmpLen<=MAXREADLENGTH ? cmpLen:MAXREADLENGTH;
		//cutSeqFromCtg(ctg2->ctgID,endOnCtg2+overlap,endOnCtg2+overlap+cmpLen-1,fastSequence,originOverlap);
		cutSeqFromRead(seqCtg2,lenCtg2,endOnCtg2+overlap,endOnCtg2+overlap+cmpLen-1,fastSequence);
		if(!bal)
			cutSeqFromRead(src_seq,rdArray[i].len,end+1,end+cmpLen,slowSequence);
		else
			cutSeqFromRead(bal_seq,rdArray[i].len,end+1,end+cmpLen,slowSequence);
		match = compareSequences(fastSequence,slowSequence, cmpLen, cmpLen);
		alignLen += cmpLen;
		matchLen += match;
		/*
		if(cmpLen>0&&match!=cmpLen+overlap){
			printSeq(stderr,fastSequence,cmpLen+overlap);
			printSeq(stderr,slowSequence,cmpLen+overlap);
			printKmer(stderr,kmerCtg2[endOnCtg2],overlap);
			fprintf(stderr,": %d(%d)\n",bal,endOnCtg2);
		}else if(cmpLen>0&&match==cmpLen+overlap)
			fprintf(stderr,"Perfect\n");
		*/
		double score = (double)matchLen/alignLen;
		if(maxScore<score){
			maxScore = score;
			//fprintf(stderr,"%4.2f (%d/%d)\n",maxScore,matchLen,alignLen);
			maxIndex = i;
		}
		SCORE[i] = score;
	}
	/*
	if(maxScore>0.0)
		fprintf(stderr,"SCORE: %4.2f\n",maxScore);
	*/
	if(maxScore>0.9){
		/*
		for(i=0;i<lenCtg1;i++)
			fprintf(stderr,"%c",int2base(seqCtg1[i]));
		fprintf(stderr,": CTG1\n");
		for(i=0;i<lenCtg2;i++)
			fprintf(stderr,"%c",int2base(seqCtg2[i]));
		fprintf(stderr,": CTG2\n");
		fprintf(stderr,"%d+%d -- %d+%d, SCORE: %4.2f\n ",offset1,offset2,cut1,cut2,maxScore);
		*/
		getSeqFromRead(rdArray[maxIndex],src_seq);
		reverseComplementSeq(src_seq, rdArray[maxIndex].len,bal_seq);

		int leftRemain = offset1-(len1-INDEX1[maxIndex]-1)>0 ? offset1-(len1-INDEX1[maxIndex]-1):0;
		int rightRemain = offset2-(overlap+INDEX2[maxIndex])>0 ? offset2-(overlap+INDEX2[maxIndex]):0;

		ctg1->gapSeqOffset = gapSeqArray->item_c;
		ctg1->gapSeqLen = END[maxIndex]-START[maxIndex]+leftRemain+rightRemain;
		if(darrayPut(gapSeqArray,ctg1->gapSeqOffset+(END[maxIndex]-START[maxIndex]+leftRemain+rightRemain)/4)){
			pt = (char *)darrayPut(gapSeqArray,ctg1->gapSeqOffset);
			for(j=0;j<leftRemain;j++){//get the left side of the gap region from search
				writeChar2tightString(seqGap[j],pt,j);
				fprintf(stderr,"%c",int2base(seqGap[j]));
			}
			for(j=START[maxIndex]+1;j<=END[maxIndex];j++){
				if(BAL[maxIndex]){
					writeChar2tightString(bal_seq[j],pt,j-START[maxIndex]-1+leftRemain);
					fprintf(stderr,"%c",int2base(bal_seq[j]));
				}else{
					writeChar2tightString(src_seq[j],pt,j-START[maxIndex]-1+leftRemain);
					fprintf(stderr,"%c",int2base(src_seq[j]));
				}
			}
			for(j=offset2-rightRemain;j<offset2;j++){   //get the right side of the gap region from search
				writeChar2tightString(seqGap[j+leftRemain],pt,j+END[maxIndex]-START[maxIndex]+leftRemain);
				fprintf(stderr,"%c",int2base(seqGap[j+leftRemain]));
			}
			
			fprintf(stderr,": GAPSEQ (%d+%d)(%d+%d)(%d+%d)(%d+%d) B %d\n",offset1,offset2,cut1,cut2,
				len1-INDEX1[maxIndex]-1,INDEX2[maxIndex],START[maxIndex],END[maxIndex],BAL[maxIndex]);
			
			ctg1->cutTail=len1-INDEX1[maxIndex]-1-offset1+cut1>cut1 ?len1-INDEX1[maxIndex]-1-offset1+cut1:cut1;
			ctg2->cutHead=overlap+INDEX2[maxIndex]-offset2+cut2>cut2 ?overlap+INDEX2[maxIndex]-offset2+cut2:cut2;
			ctg2->scaftig_start = 0;
			ret = 1;
		}
	}
	free((void*)START);
	free((void*)END);
	free((void*)INDEX1);
	free((void*)INDEX2);
	free((void*)SCORE);
	free((void*)BAL);

	free((void*)src_seq);
	free((void*)bal_seq);
	return ret;
}


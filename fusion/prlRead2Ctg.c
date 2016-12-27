#include "stdinc.h" 
#include "newhash.h"
#include "extfunc.h" 
#include "extvab.h" 

static ubyte binLight[16384];
static ubyte probableMatrix[4][4] = {
// A  C  T  G								 7 3 2 1
	{7, 2, 1, 3},				//	A->A G C T
	{2, 7, 3, 1},				//	C->C T A G
	{1, 3, 7, 2},				//	T->T C G A
	{3, 1, 2, 7}				//	G->G A T C	
};

static ubyte2 doubleBitMasker[7] = {
	0x3,								//000000 00000011
	0xC,								//000000 00001100
	0x30,								//000000 00110000
	0xC0,								//000000 11000000
	0x300,							//000011 00000000
	0xC00,							//001100 00000000
	0x3000							//110000 00000000
};

static boolean staticFlag=1;

static long long readsInGap=0;

static int buffer_size=10000000;
static long long readCounter;
static long long mapCounter;
long long single_count;
long long single_map;
static int ALIGNLEN=0;
//buffer related varibles for chop kmer
static int read_c;
static char **rcSeq;
static char **seqBuffer;
static int *lenBuffer;
static unsigned int *ctgIdArray;
static int *posArray;
static char *orienArray;
static char *footprint;   // flag indicates whether the read shoulld leave markers on contigs

// kmer related variables
static int kmer_c;
static Kmer *kmerBuffer,*hashBanBuffer;
static kmer_t **nodeBuffer;
static boolean *smallerBuffer;
static unsigned int *indexArray;

static int *deletion;

static void parse1read(int t,int threadID);
static void threadRoutine(void *thrdID);
static void searchKmer(int t,KmerSet *kset);
static void chopKmer4read(int t,int threadID);
static void thread_wait(pthread_t *threads);

static void creatThrds(pthread_t *threads,PARAMETER *paras)
{
	unsigned char i;
	int temp;
	
	for(i=0;i<thrd_num;i++){
		//printf("to create %dth thread\n",(*(char *)&(threadID[i])));
		if((temp = pthread_create(&threads[i], NULL, (void *)threadRoutine, &(paras[i]))) != 0){ 
			printf("create threads failed\n");
			exit(1);
		}
	}
	//printf("%d thread created\n",thrd_num);
}

static void threadRoutine(void *para)
{
	PARAMETER *prm;
	int i,t;
	unsigned char id;
	
	prm = (PARAMETER *)para;
	id = prm->threadID;
	//printf("%dth thread with task %d, hash_table %p\n",id,prm.task,prm.hash_table);
	while(1){
		if(*(prm->selfSignal)==1){
			for(i=0;i<kmer_c;i++)
			{
				//if((hashBanBuffer[i]&taskMask)!=prm.threadID)
				if((hashBanBuffer[i]%thrd_num)!=id)
					continue;
				searchKmer(i,KmerSets[id]);
			}
			*(prm->selfSignal) = 0;
		}else if(*(prm->selfSignal)==2){
			for(i=0;i<read_c;i++){
				if(i%thrd_num!=id)
					continue;
				chopKmer4read(i,id+1);
			}
			*(prm->selfSignal) = 0;
		}else if(*(prm->selfSignal)==3){
		// parse reads
			for(t=0;t<read_c;t++){
				if(t%thrd_num!=id)
					continue;
				parse1read(t,id+1);
			}
			*(prm->selfSignal) = 0;
		}else if(*(prm->selfSignal)==5){
			*(prm->selfSignal) = 0;
			break;
		}

		usleep(1);
	}
}
/*
static void chopReads()
{
	int i;
	for(i=0;i<read_c;i++){
		chopKmer4read(i,0);
	}
}
*/
static void chopKmer4read(int t,int threadID)
{
	int len_seq=lenBuffer[t];

	if(len_seq<overlaplen+1)
		return;
	
	char *src_seq=seqBuffer[t];
	char *bal_seq=rcSeq[threadID];
	int j,bal_j;
	Kmer hash_ban,bal_hash_ban;
	Kmer word,bal_word;
	int index;
	
	word = 0;
	for (index = 0;index<overlaplen;index++){
		word <<= 2;
		word += src_seq[index];
	}
	reverseComplementSeq(src_seq, len_seq,bal_seq);

		// complementary node
	bal_word = reverseComplement(word,overlaplen);
	bal_j = len_seq-0-overlaplen;  //  0;     
	index = indexArray[t];
	if(word<bal_word){
		hash_ban = hash_kmer(word);
		kmerBuffer[index] = word;
		smallerBuffer[index] = 1;
		hashBanBuffer[index++] = hash_ban;
	}else{
		bal_hash_ban = hash_kmer(bal_word);
		kmerBuffer[index] = bal_word;
		smallerBuffer[index] = 0;
		hashBanBuffer[index++] = bal_hash_ban;
	}

		//printf("%dth: %p with %p\n",kmer_c-1,bal_word,bal_hash_ban);
	for(j = 1; j <= len_seq - overlaplen; j ++)	{
		word = nextKmer(word,src_seq[j-1+overlaplen]);
		bal_j = len_seq-j-overlaplen; //  j; 
		bal_word = prevKmer(bal_word,bal_seq[bal_j]);
		if(word<bal_word){
			hash_ban = hash_kmer(word);
			kmerBuffer[index] = word;
			smallerBuffer[index] = 1;
			hashBanBuffer[index++] = hash_ban;
				//printf("%dth: %p with %p\n",kmer_c-1,word,hashBanBuffer[kmer_c-1]);
		}else{
				// complementary node
			bal_hash_ban = hash_kmer(bal_word); 
			kmerBuffer[index] = bal_word;
			smallerBuffer[index] = 0;
			hashBanBuffer[index++] = bal_hash_ban;
				//printf("%dth: %p with %p\n",kmer_c-1,bal_word,hashBanBuffer[kmer_c-1]);
		}
	}
}

//splay for one kmer in buffer and save the node to nodeBuffer
static void searchKmer(int t,KmerSet *kset)
{
	kmer_t *node;
	
	boolean found = search_kmerset(kset, kmerBuffer[t], &node);

	++kset->searchCnt;
	
	if(found)
	{
		++kset->foundCnt;
		if(!node->deleted)
			nodeBuffer[t] = node;
		else
		{
			++kset->delCnt;
			nodeBuffer[t] = NULL;
		}
	}
	else
	{
		++kset->searchSpcSeedCnt;

		boolean spcFlag;
		Kmer buff_kmer, spc_kmer;
		ubyte2 spc_bases;
		spcKmer *rs;
		spcBase *tmpBase;
		
		buff_kmer=kmerBuffer[t];
		spc_kmer = ((buff_kmer>>14)&0xFFFFFFF00) | ((buff_kmer>>12)&0xC0) | ((buff_kmer>>10)&0x3C) | ((buff_kmer>>6)&0x3);
		spc_bases = ((buff_kmer>>8)&0x3000) | ((buff_kmer>>6)&0xC00) | ((buff_kmer>>2)&0x3C0) | (buff_kmer&0x3F);
		
		spcFlag = search_spckmerset(spcSet, spc_kmer, &rs);
		
		if(spcFlag)
		{
			++kset->getSpcSeedCnt;
			
			int i=0,j=0,getFlag=-1;
			int mismatch=0;
			ubyte2 tmp,mostLastBase;	//loci flags
			ubyte2 bestSpcBases;			//best spaced bases
			int min_mis=31;
			ubyte2 tmpSpcBase;
							
			tmpBase=rs->start;
			
			//fprintf(stderr,"search %llu\tspc_kmer %u\tspc_bases %u\n", kmerBuffer[t], spc_kmer, spc_bases);
			
			while(tmpBase != NULL)
			{
				tmpSpcBase = tmpBase->spaced_bases;
				tmp = ((spc_bases ^ tmpSpcBase) & 0x5555) | (((spc_bases ^ tmpSpcBase) & 0xAAAA)>>1);
				mismatch=binLight[tmp];
				
				if(mismatch < min_mis)				//get the minimal mismatch spaced_bases
				{
					min_mis = mismatch;
					mostLastBase = tmp;
					bestSpcBases = tmpSpcBase;
					node = tmpBase->large_kmer;
					getFlag=0;
				}
				else if(mismatch == min_mis)	//if same amount of mismatch, choose the most right mismatch pos
				{
					if(tmp<mostLastBase)
					{
						mostLastBase = tmp;
						bestSpcBases = tmpSpcBase;
						node = tmpBase->large_kmer;
						getFlag=1;
					}
					else if(tmp == mostLastBase)	//if same mismatch pos, choose the most probable one[see probableMatrix]
					{
/*
static ubyte probableMatrix[4][4] = {
//A  C  T  G							 7 3 2 1
7, 2, 1, 3,				//	A->A G C T
2, 7, 3, 1,				//	C->C T A G
1, 3, 7, 2,				//	T->T C G A
3, 1, 2, 7				//	G->G A T C	
};
*/
						getFlag=2;
						ubyte2 readBases = spc_bases, loopBases = tmpSpcBase, bestBases = bestSpcBases, mismatchFlag = tmp;
						for(j=0;j<7;j++)
						{
							if((mismatchFlag & 0x3) > 0)
							{
								if(probableMatrix[(readBases & 0x3)][(loopBases & 0x3)] > probableMatrix[(readBases & 0x3)][(bestBases & 0x3)])
									//check each 2 bits(1 base) if mismatch
								{
									mostLastBase = tmp;
									bestSpcBases = tmpSpcBase;
									node = tmpBase->large_kmer;
									break;
								}
								else if((probableMatrix[(readBases & 0x3)][(loopBases & 0x3)] < probableMatrix[(readBases & 0x3)][(bestBases & 0x3)]))
									break;
							}
							mismatchFlag>>=2;
							readBases>>=2;
							loopBases>>=2;
							bestBases>>=2;
						}
					}
				}
				
				tmpBase = tmpBase->next;
			}

			if(getFlag<0)
			{
				fprintf(stderr,"getFlag error at %llu",kmerBuffer[t]);
				exit(-1);
			}
			++kset->levelGet[getFlag];
			nodeBuffer[t] = node;
		}
		else
			nodeBuffer[t] = NULL;
	}
}

static void parse1read(int t,int threadID)
{
	unsigned int j,i,s;
	unsigned int contigID;
	int counter2=0,counter;
	unsigned int ctgLen,pos;
	kmer_t *node;
	boolean isSmaller;
	int flag,maxOcc=0;
	kmer_t *maxNode=NULL;
	int alldgnLen = lenBuffer[t] > ALIGNLEN ? ALIGNLEN:lenBuffer[t];
	int multi = alldgnLen-overlaplen+1 < 5 ? 5:alldgnLen-overlaplen+1;
	unsigned int start,finish;

	footprint[t] = 0;

	start = indexArray[t];
	finish = indexArray[t+1];
	if(finish==start){  //too short
		ctgIdArray[t] = 0;
		deletion[threadID]++;
		return;
	}
	for(j=start;j<finish;j++)
		if(nodeBuffer[j])
			counter2++;
	if(counter2<2)
		deletion[threadID]++;
	counter = counter2 = 0;
	for(j=start;j<finish;j++){
		node = nodeBuffer[j];	
		if(!node)  //same as previous
			continue;
		flag = 1;
		for(s=j+1;s<finish;s++){
			if(!nodeBuffer[s])
				continue;
			if(nodeBuffer[s]->l_links==node->l_links){
				flag++;
				nodeBuffer[s] = NULL;
			}
		}
		if(flag>=2)
			counter2++; //a loose alignment 
		if(flag>=multi)
			counter++;
		else
			continue;
		if(flag>maxOcc){
			pos = j;
			maxOcc = flag;
			maxNode = node;
		}
	}
	if(!counter){  //no match
		ctgIdArray[t] = 0;
		return;
	}
	if(counter2>1)
		footprint[t] = 1;  //aligned to multi contigs

	j = pos;
	i = pos - start + 1;
	node = nodeBuffer[j];	
	isSmaller = smallerBuffer[j];
	contigID = node->l_links;
	ctgLen = contig_array[contigID].length;
	pos = node->r_links;
	if(node->twin==isSmaller){
		orienArray[t] = '-';
		ctgIdArray[t] = getTwinCtg(contigID);
		posArray[t] = ctgLen - pos -overlaplen -i + 1;
	}else{
		orienArray[t] = '+';
		ctgIdArray[t] = contigID;
		posArray[t] = pos - i + 1;
	}
	
}

static void sendWorkSignal(unsigned char SIG,unsigned char *thrdSignals)
{
	int t;

	for(t=0;t<thrd_num;t++)
		thrdSignals[t+1] = SIG;
	while(1){
		usleep(10);
				
		for(t=0;t<thrd_num;t++)
			if(thrdSignals[t+1])
				break;
		if(t==thrd_num)
			break;
	}
}

static void locate1read(int t)
{
	int i,j,start,finish;
	kmer_t *node;
	unsigned int contigID;
	int pos,ctgLen;
	boolean isSmaller;

	start = indexArray[t];
	finish = indexArray[t+1];
		
	for(j=start;j<finish;j++){
		node = nodeBuffer[j];	
		if(!node)  //same as previous
			continue;
		i = j - start + 1;
		isSmaller = smallerBuffer[j];
		contigID = node->l_links;
		ctgLen = contig_array[contigID].length;
		pos = node->r_links;
		if(node->twin==isSmaller){
			ctgIdArray[t] = getTwinCtg(contigID);
			posArray[t] = ctgLen - pos -overlaplen -i + 1;
		}else{
			ctgIdArray[t] = contigID;
			posArray[t] = pos - i + 1;
		}
	}
	
}

static void output1read(int t, FILE *outfp)
{
	int len = lenBuffer[t];
	int index;
	readsInGap++;
/*
	if(ctgIdArray[t]==735||ctgIdArray[t]==getTwinCtg(735)){
		printf("%d\t%d\t%d\t",t+1,ctgIdArray[t],posArray[t]);
		int j;
		for(j=0;j<len;j++)
			printf("%c",int2base((int)seqBuffer[t][j]));
		printf("\n");
	}
*/
	for(index=0;index<len;index++)
		writeChar2tightString(seqBuffer[t][index],rcSeq[1],index);
	fwrite(&len,sizeof(int),1,outfp);
	fwrite(&ctgIdArray[t],sizeof(int),1,outfp);
	fwrite(&posArray[t],sizeof(int),1,outfp);
	fwrite(rcSeq[1],sizeof(char),len/4+1,outfp);

}

static void getReadIngap(int t,int insSize,FILE *outfp,boolean readOne)
{
	int read1,read2;
	if(readOne){
		read1 = t;
		read2 = t + 1;
		ctgIdArray[read1] = ctgIdArray[read2];
		posArray[read1] = posArray[read2] + insSize - lenBuffer[read1];  //   --> R2       <-- R1
		output1read(read1,outfp);
	}else{
		read2 = t;
		read1 = t - 1;
		ctgIdArray[read2] = ctgIdArray[read1];
		posArray[read2] = posArray[read1] + insSize - lenBuffer[read2];  // --> R1     <-- R2
		output1read(read2,outfp);
	}
}

static void recordLongRead(FILE *outfp)
{
	int t;

	for(t=0;t<read_c;t++){
		readCounter++;
		if(footprint[t])
			output1read(t,outfp);
	}
}

static void recordAlldgn(FILE *outfp,int insSize,FILE *outfp2)
{
	int t,ctgId;
	boolean rd1gap,rd2gap;

	for(t=0;t<read_c;t++){
		readCounter++;
		single_count++;
		rd1gap = rd2gap = 0;
		ctgId = ctgIdArray[t];
		if(outfp2&&t%2==1){ //make sure this is read2 in a pair
			if(ctgIdArray[t]<1&&ctgIdArray[t-1]>0){
				getReadIngap(t,insSize,outfp2,0);  //read 2 in gap
				rd2gap = 1;
			}
			else if(ctgIdArray[t]>0&&ctgIdArray[t-1]<1){
				getReadIngap(t-1,insSize,outfp2,1);  //read 1 in gap
				rd1gap = 1;
			}
		}
		if(ctgId<1)
			continue;
		mapCounter++;
		single_map++;
		fprintf(outfp,"%lld\t%u\t%d\t%c\n",readCounter,
						ctgIdArray[t],posArray[t],orienArray[t]);
		if(t%2==0)
			continue;
		if(outfp2&&footprint[t-1]&&!rd1gap)
			output1read(t-1,outfp2);
		if(outfp2&&footprint[t]&&!rd2gap)
			output1read(t,outfp2);
		
	}
}

//load contig index and length
void basicContigInfo(char *infile)
{
	char name[256],lldne[1024];
	FILE *fp;
	int length,bal_ed,num_all,num_long,index;
	
	sprintf(name,"%s.ContigIndex",infile);
	fp = ckopen(name,"r");

	fgets(lldne,sizeof(lldne),fp);
	sscanf(lldne+8,"%d %d",&num_all,&num_long);
	//printf("%d edges in graph\n",num_all);
	num_ctg = num_all;
	contig_array = (CONTIG *)ckalloc((num_all+1)*sizeof(CONTIG));

	fgets(lldne,sizeof(lldne),fp);
	num_long = 0;
	while(fgets(lldne,sizeof(lldne),fp)!=NULL){
		sscanf(lldne,"%d %d %d",&index,&length,&bal_ed);
		
		contig_array[++num_long].length = length;
		contig_array[num_long].bal_edge = bal_ed+1;
		if(index!=num_long)
			printf("basicContigInfo: %d vs %d\n",index,num_long);
		if(bal_ed==0)
			continue;
		contig_array[++num_long].length = length;
		contig_array[num_long].bal_edge = -bal_ed+1;
			
	}

	fclose(fp);
}

void prlRead2Ctg(char *libfile,char *outfile)
{
	long long i;
	char *src_name,*next_name,name[256];
	FILE *fo,*outfp2=NULL;
	int maxReadNum,libNo,prevLibNo,insSize;
	boolean flag,pairs=1;
	pthread_t threads[thrd_num];
	unsigned char thrdSignal[thrd_num+1];
	PARAMETER paras[thrd_num];

	maxReadLen = 0;
	maxNameLen = 256;
	scan_libInfo(libfile);
	alloc_pe_mem(num_libs);
	if(!maxReadLen)
		maxReadLen = 100;
	//printf("In file: %s, max seq len %d, max name len %d\n\n",
		//libfile,maxReadLen,maxNameLen);
	if(maxReadLen>maxReadLen4all)
		maxReadLen4all = maxReadLen;

//////////////////////////////////////////// spcSet
	fflush(stdout);

	ubyte2 spc_i,spc_j;
	for(spc_i=0;spc_i<16384;spc_i++)
	{
		binLight[spc_i]=0;
		for(spc_j=spc_i;spc_j;spc_j=spc_j&(spc_j-1))
			++binLight[spc_i];
	}

	spcSet = init_spckmerset(KmerSets[thrd_num-1]->size*thrd_num, 0.77f);

	
	for(i=0;i<thrd_num;i++)
	{
		buildSpcKmerSet(KmerSets[i], spcSet);
		mvnv(0,"%lldth spaced bases set build complete.",i);
	}
//////////////////////////////////////////// END spcSet

	src_name = (char *)ckalloc((maxNameLen+1)*sizeof(char));
	next_name = (char *)ckalloc((maxNameLen+1)*sizeof(char));

	kmerBuffer = (Kmer *)ckalloc(buffer_size*sizeof(Kmer));
	hashBanBuffer = (Kmer *)ckalloc(buffer_size*sizeof(Kmer));
	nodeBuffer = (kmer_t **)ckalloc(buffer_size*sizeof(kmer_t *));
	smallerBuffer = (boolean *)ckalloc(buffer_size*sizeof(boolean));

	maxReadNum = buffer_size/(maxReadLen-overlaplen+1);
	maxReadNum = maxReadNum%2==0 ? maxReadNum:maxReadNum-1;   //make sure paired reads are processed at the same batch
	seqBuffer = (char **)ckalloc(maxReadNum*sizeof(char*));
	lenBuffer = (int *)ckalloc(maxReadNum*sizeof(int));
	indexArray = (unsigned int *)ckalloc((maxReadNum+1)*sizeof(unsigned int));
	ctgIdArray = (unsigned int *)ckalloc((maxReadNum+1)*sizeof(unsigned int));
	posArray = (int *)ckalloc((maxReadNum+1)*sizeof(int));
	orienArray = (char *)ckalloc((maxReadNum+1)*sizeof(char));
	footprint = (char *)ckalloc((maxReadNum+1)*sizeof(char));

	for(i=0;i<maxReadNum;i++)
		seqBuffer[i] = (char *)ckalloc(maxReadLen*sizeof(char));

	rcSeq = (char **)ckalloc((thrd_num+1)*sizeof(char*));
	deletion = (int *)ckalloc((thrd_num+1)*sizeof(int));
	thrdSignal[0] = 0;
	deletion[0] = 0;
	if(1){
		for(i=0;i<thrd_num;i++){
			rcSeq[i+1] = (char *)ckalloc(maxReadLen*sizeof(char)); 
			deletion[i+1] = 0;
			thrdSignal[i+1] = 0;
			paras[i].threadID = i;
			paras[i].mainSignal = &thrdSignal[0];
			paras[i].selfSignal = &thrdSignal[i+1];
		}
		creatThrds(threads,paras);
	}
	
	if(!contig_array)
		basicContigInfo(outfile);
	
	sprintf(name,"%s.readInGap",outfile);
	outfp2 = ckopen(name, "wb");

	sprintf(name,"%s.readOnContig",outfile);
	fo = ckopen(name,"w");
	fprintf(fo,"read\tcontig\tpos\n");
	readCounter = mapCounter = 0;
	single_map = single_count = 0;
	kmer_c = n_solexa = read_c = i = libNo = readNumBack = gradsCounter = readsInGap = 0;
	prevLibNo = -1;
	
	//mvnv(0,"Start loading reads.");
	
	while((flag=read1seqInLib(seqBuffer[read_c],next_name,&(lenBuffer[read_c]),&libNo,pairs,0))!=0)	{
		if(libNo!=prevLibNo){
			prevLibNo = libNo;
			insSize = lib_array[libNo].avg_ins;
			ALIGNLEN = lib_array[libNo].map_len;
			
			if(insSize>1000)
				ALIGNLEN = ALIGNLEN < 35 ? 35: ALIGNLEN;
			else
				ALIGNLEN = ALIGNLEN < 32 ? 32: ALIGNLEN;
			//printf("current insert size %d, map_len %d\n",insSize,ALIGNLEN);
			
		}
	
		if(insSize>1000)  
			ALIGNLEN = ALIGNLEN < (lenBuffer[read_c]/2+1) ? (lenBuffer[read_c]/2+1):ALIGNLEN;
		
//		if((++i)%100000000==0)
//			printf("[%s]%lld reads processed.\n",__FUNCTION__,i);
		indexArray[read_c] = kmer_c;
		if(lenBuffer[read_c] >= overlaplen+1)     
			kmer_c += lenBuffer[read_c] - overlaplen + 1;
		read_c++;
		if(read_c==maxReadNum){
			//mvnv(0,"Start processing reads.");
			
			indexArray[read_c] = kmer_c;
			
			sendWorkSignal(2,thrdSignal);
			//mvnv(0,"chop finished one buffer.");
			sendWorkSignal(1,thrdSignal);
			//mvnv(0,"search finished one buffer.");
			sendWorkSignal(3,thrdSignal);
			//mvnv(0,"parse finished one buffer.");

			recordAlldgn(fo,insSize,outfp2);
			kmer_c = 0;
			read_c = 0;
		}
	}

	if(read_c){
		indexArray[read_c] = kmer_c;
		sendWorkSignal(2,thrdSignal);
		sendWorkSignal(1,thrdSignal);
		sendWorkSignal(3,thrdSignal);
		recordAlldgn(fo,insSize,outfp2);
		//printf("Output %lld out of %lld (%.1f)%% reads in gaps\n",readsInGap,readCounter,
			//		(float)readsInGap/readCounter*100);
	}
	if(readCounter)
		printf("[%s]total %llu reads , map-rate (%.1f)%%\n",__FUNCTION__,
			readCounter,(float)mapCounter/readCounter*100);
	sendWorkSignal(5,thrdSignal);
	
	thread_wait(threads);
	fclose(fo);

	sprintf(name,"%s.peGrads",outfile);	
	fo = ckopen(name,"w");
	fprintf(fo,"grads&num: %d\t%lld\t%d\n",gradsCounter,n_solexa,maxReadLen4all);
	if(pairs){ 
		if(gradsCounter)
			;
			//printf("%d pe insert size, the largest boundary is %lld\n\n",
				//gradsCounter,pes[gradsCounter-1].PE_bound);
		else
			printf("no paired reads found\n");
		for(i=0;i<gradsCounter;i++)
			fprintf(fo,"%d\t%lld\t%d\t%d\n",pes[i].insertS,pes[i].PE_bound,pes[i].rank,pes[i].pair_num_cut);
	}
	fclose(fo);

	fclose(outfp2);
	
	free_pe_mem();
	free_libs();
	if(1){ // multi-threads
		for(i=0;i<thrd_num;i++){
			deletion[0] += deletion[i+1];
			free((void *)rcSeq[i+1]);
		}
	}
	//printf("%d reads deleted\n",deletion[0]);
	
	ubyte8 searchCntTot=0, foundCntTot=0, delCntTot=0,
					searchSpcSeedCntTot=0,getSpcSeedCntTot=0,levelGet1=0,levelGet2=0,levelGet3=0;
	for(i=0;i<thrd_num;i++){
		searchCntTot += KmerSets[i]->searchCnt;
		foundCntTot += KmerSets[i]->foundCnt;
		delCntTot += KmerSets[i]->delCnt;
		searchSpcSeedCntTot += KmerSets[i]->searchSpcSeedCnt;
		getSpcSeedCntTot += KmerSets[i]->getSpcSeedCnt;
		levelGet1 += KmerSets[i]->levelGet[0];
		levelGet2 += KmerSets[i]->levelGet[1];
		levelGet3 += KmerSets[i]->levelGet[2];		
	}
	fprintf(stderr,"SEARCH: Search %llu, get %llu, deleted %llu\n", 
									searchCntTot, foundCntTot, delCntTot);	
	fprintf(stderr,"SPACED SEED: Search %llu, get %llu, LVnum %llu, LVpos %llu, LVpro %llu\n", 
									searchSpcSeedCntTot, getSpcSeedCntTot, levelGet1, levelGet2, levelGet3);
	
	free((void *)rcSeq);
	free((void *)deletion);
	for(i=0;i<maxReadNum;i++)
		free((void *)seqBuffer[i]);
	free((void *)seqBuffer);
	free((void *)lenBuffer);
	free((void *)indexArray);
	
	free((void *)kmerBuffer);
	free((void *)smallerBuffer);
	free((void *)hashBanBuffer);
	free((void *)nodeBuffer); 
	free((void *)ctgIdArray);
	free((void *)posArray);
	free((void *)orienArray);
	free((void *)footprint);

	free((void *)src_name);
	free((void *)next_name);

	if(contig_array){
		free((void *)contig_array);
		contig_array = NULL;
	}
}

static void thread_wait(pthread_t *threads)
{
	int i;
	for(i=0;i<thrd_num;i++)
		if(threads[i]!=0)
			pthread_join(threads[i],NULL);
	
}

/********************* map long reads for gap filling ************************/
void prlLongRead2Ctg(char *libfile,char *outfile)
{
	long long i;
	char *src_name,*next_name,name[256];
	FILE *outfp2;
	int maxReadNum,libNo,prevLibNo;
	boolean flag,pairs=0;
	pthread_t threads[thrd_num];
	unsigned char thrdSignal[thrd_num+1];
	PARAMETER paras[thrd_num];

	maxReadLen = 0;
	maxNameLen = 256;
	scan_libInfo(libfile);
	if(!maxReadLen)
		maxReadLen = 100;
	int longReadLen = getMaxLongReadLen(num_libs);
	if(longReadLen<1)  // no long reads
		return;
	maxReadLen4all = maxReadLen<longReadLen ? longReadLen:maxReadLen;
	//printf("In file: %s, long read len %d, max name len %d\n\n",
		//libfile,longReadLen,maxNameLen);
	maxReadLen = longReadLen;

	src_name = (char *)ckalloc((maxNameLen+1)*sizeof(char));
	next_name = (char *)ckalloc((maxNameLen+1)*sizeof(char));

	kmerBuffer = (Kmer *)ckalloc(buffer_size*sizeof(Kmer));
	hashBanBuffer = (Kmer *)ckalloc(buffer_size*sizeof(Kmer));
	nodeBuffer = (kmer_t **)ckalloc(buffer_size*sizeof(kmer_t *));
	smallerBuffer = (boolean *)ckalloc(buffer_size*sizeof(boolean));

	maxReadNum = buffer_size/(maxReadLen-overlaplen+1);
	maxReadNum = maxReadNum%2==0 ? maxReadNum:maxReadNum-1;   //make sure paired reads are processed at the same batch
	seqBuffer = (char **)ckalloc(maxReadNum*sizeof(char*));
	lenBuffer = (int *)ckalloc(maxReadNum*sizeof(int));
	indexArray = (unsigned int *)ckalloc((maxReadNum+1)*sizeof(unsigned int));
	ctgIdArray = (unsigned int *)ckalloc((maxReadNum+1)*sizeof(unsigned int));
	posArray = (int *)ckalloc((maxReadNum+1)*sizeof(int));
	orienArray = (char *)ckalloc((maxReadNum+1)*sizeof(char));
	footprint = (char *)ckalloc((maxReadNum+1)*sizeof(char));

	for(i=0;i<maxReadNum;i++)
		seqBuffer[i] = (char *)ckalloc(maxReadLen*sizeof(char));

	rcSeq = (char **)ckalloc((thrd_num+1)*sizeof(char*));
	deletion = (int *)ckalloc((thrd_num+1)*sizeof(int));
	thrdSignal[0] = 0;
	deletion[0] = 0;
	if(1){
		for(i=0;i<thrd_num;i++){
			rcSeq[i+1] = (char *)ckalloc(maxReadLen*sizeof(char)); 
			deletion[i+1] = 0;
			thrdSignal[i+1] = 0;
			paras[i].threadID = i;
			paras[i].mainSignal = &thrdSignal[0];
			paras[i].selfSignal = &thrdSignal[i+1];
		}
		creatThrds(threads,paras);
	}
	
	if(!contig_array)
		basicContigInfo(outfile);
	
	sprintf(name,"%s.longReadInGap",outfile);
	outfp2 = ckopen(name, "wb");

	readCounter = 0;
	kmer_c = n_solexa = read_c = i = libNo;
	prevLibNo = -1;
	while((flag=read1seqInLib(seqBuffer[read_c],next_name,&(lenBuffer[read_c]),&libNo,pairs,4))!=0)	{
		if(libNo!=prevLibNo){
			prevLibNo = libNo;
			ALIGNLEN = lib_array[libNo].map_len;
			ALIGNLEN = ALIGNLEN < 35 ? 35: ALIGNLEN;
			//printf("Map_len %d\n",ALIGNLEN);
		}
	
//		if((++i)%100000000==0)
	//		printf("%lld reads processed.\n",i);
		indexArray[read_c] = kmer_c;
		if(lenBuffer[read_c] >= overlaplen+1)     
			kmer_c += lenBuffer[read_c] - overlaplen + 1;
		read_c++;
		if(read_c==maxReadNum){
			indexArray[read_c] = kmer_c;
			
			sendWorkSignal(2,thrdSignal);
			sendWorkSignal(1,thrdSignal);
			sendWorkSignal(3,thrdSignal);

			recordLongRead(outfp2);
			kmer_c = 0;
			read_c = 0;
		}
	}

	if(read_c){
		indexArray[read_c] = kmer_c;
		sendWorkSignal(2,thrdSignal);
		sendWorkSignal(1,thrdSignal);
		sendWorkSignal(3,thrdSignal);
		recordLongRead(outfp2);
		//printf("Output %lld out of %lld (%.1f)%% reads in gaps\n",readsInGap,readCounter,
			//		(float)readsInGap/readCounter*100);
	}

	sendWorkSignal(5,thrdSignal);
	
	thread_wait(threads);

	fclose(outfp2);
	
	free_libs();
	if(1){ // multi-threads
		for(i=0;i<thrd_num;i++){
			deletion[0] += deletion[i+1];
			free((void *)rcSeq[i+1]);
		}
	}
	//printf("%d reads deleted\n",deletion[0]);
	
	free((void *)rcSeq);
	free((void *)deletion);
	for(i=0;i<maxReadNum;i++)
		free((void *)seqBuffer[i]);
	free((void *)seqBuffer);
	free((void *)lenBuffer);
	free((void *)indexArray);
	
	free((void *)kmerBuffer);
	free((void *)smallerBuffer);
	free((void *)hashBanBuffer);
	free((void *)nodeBuffer); 
	free((void *)ctgIdArray);
	free((void *)posArray);
	free((void *)orienArray);
	free((void *)footprint);

	free((void *)src_name);
	free((void *)next_name);

}


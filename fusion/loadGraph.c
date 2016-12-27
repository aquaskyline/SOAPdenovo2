#include "stdinc.h" 
#include "newhash.h"
#include "extfunc.h" 
#include "extvab.h" 

#define preARCBLOCKSIZE 100000

static void loadArcs(char *graphfile);
static void loadContig(char *graphfile);

void loadUpdatedVertex(char *graphfile)
{
	char name[256],line[256];
	FILE *fp;
	Kmer word,bal_word;
	int num_kmer,i;
	char ch;

	sprintf(name,"%s.updated.vertex",graphfile);
	fp = ckopen(name,"r");

	while(fgets(line,sizeof(line),fp)!=NULL){
		if(line[0] == 'V'){
			sscanf(line+6, "%d %c %d",&num_kmer,&ch,&overlaplen);
			printf("there're %d kmers in vertex file\n",num_kmer);
			//printf("total %d kmer in all contigs.\n",num_kmer);
			break;
		}
	}
	
	vt_array = (VERTEX *)ckalloc((2*num_kmer)*sizeof(VERTEX));
	for(i=0;i<num_kmer;i++){
		fscanf(fp,"%llx ",&word);
		vt_array[i].kmer = word;
	}
	fclose(fp);
	
	for(i=0;i<num_kmer;i++){
		bal_word = reverseComplement(vt_array[i].kmer,overlaplen);
		vt_array[i+num_kmer].kmer = bal_word;
	}
	num_vt = num_kmer;
}
int cmp_int(const void *a,const void *b)
{
        int *A,*B;
        A = (int *)a;
        B = (int *)b;

        if(*A>*B)
                return 1;
        else if(*A==*B)
                return 0;
        else
                return -1;
}
int uniqueLenSearch(unsigned int *len_array,unsigned int *flag_array,int num,unsigned int target)
{
	int mid,low,high;
	low = 1;
	high = num;

	while(low<=high){
		mid = (low+high)/2;
		if(len_array[mid]==target)
			break;
		else if(target>len_array[mid])
			low = mid+1;
		else
			high = mid-1;
	}
	if(low>high)
		return -1;
	//locate the first same length unflaged 
	return flag_array[mid]++;
	
}

int lengthSearch(unsigned int *len_array,unsigned int *flag_array,int num,unsigned int target)
{
	int mid,low,high,i;
	low = 1;
	high = num;

	while(low<=high){
		mid = (low+high)/2;
		if(len_array[mid]==target)
			break;
		else if(target>len_array[mid])
			low = mid+1;
		else
			high = mid-1;
	}
	if(low>high)
		return -1;
	//locate the first same length unflaged 
	if(!flag_array[mid]){
		for(i=mid-1;i>0;i--){
			if(len_array[i]!=len_array[mid]||flag_array[i])
				break;
		}
		flag_array[i+1] = 1;
		return i+1;
	}else{
		for(i=mid+1;i<=num;i++){
			if(!flag_array[i])
				break;
		}
		flag_array[i] = 1;
		return i;
	}

}

void quick_sort_int(unsigned int *length_array, int low, int high)
{
	int i, j;
	Kmer pivot;
	if (low < high)
	{
		pivot=length_array[low];
		i=low;
		j=high;

		while(i<j)
		{
			while (i<j && length_array[j]>=pivot)
			j--;
			if(i<j)
				length_array[i++]=length_array[j]; 

			while (i<j && length_array[i]<=pivot)
			i++;
			if(i<j)
			length_array[j--]=length_array[i]; 
		}

		length_array[i]=pivot;

		quick_sort_int(length_array,low,i-1);
		quick_sort_int(length_array,i+1,high);
	}
} 

void loadUpdatedEdges(char *graphfile)
{
	char c,name[256],line[1024];
	int bal_ed,cvg;
	FILE *fp,*out_fp;
	unsigned long long from_kmer,to_kmer;
	unsigned int num_ctgge,length,index=0,num_kmer;
	unsigned int i=0,j;
	unsigned int newIndex;
	unsigned int *length_array,*flag_array,diff_len;
	char *outfile = graphfile;
	long long cvgSum=0;
	long long counter=0;

	//get overlaplen from *.preGraphBasic
	/*sprintf(name, "%s.preGraphBasic", graphfile);
	fp = ckopen(name, "r");

	while(fgets(line,sizeof(line),fp)!=NULL){
		if(line[0] == 'V'){
			sscanf(line+6, "%d %c %d",&num_kmer,&c,&overlaplen);
			//printf("K = %d\n",overlaplen);
			break;
		}
	}*/
	if(ctg_short==0)
		ctg_short = overlaplen + 2;

	//fclose(fp);

	sprintf(name,"%s.updated.edge",graphfile);
	fp = ckopen(name,"r");
	sprintf(name,"%s.newContigIndex",outfile);
	out_fp = ckopen(name,"w");

	while(fgets(line,sizeof(line),fp)!=NULL){
		if(line[0] == 'E'){
			sscanf(line+5, "%d",&num_ctgge);
			//printf("there're %d edge in edge file\n",num_ctgge);
			//printf("total %d contigs\n",num_ctgge);
			break;
		}
	}

	index_array = (unsigned int *)ckalloc((num_ctgge+1)*sizeof(unsigned int));
	length_array = (unsigned int *)ckalloc((num_ctgge+1)*sizeof(unsigned int));
	flag_array = (unsigned int *)ckalloc((num_ctgge+1)*sizeof(unsigned int));
	while(fgets(line,sizeof(line),fp)!=NULL){
		if(line[0]=='>'){
			sscanf(line+7,"%d",&length);
			index_array[++index] = length;
			length_array[++i] = length;
		}
	}
	num_ctg = index;
	orig2new = 1;
	//quick_sort_int(length_array,1,num_ctg);
	qsort(&(length_array[1]),num_ctg,sizeof(length_array[0]),cmp_int);
	//extract unique length
	diff_len = 0;
	for(i=1;i<=num_ctg;i++){
		for(j=i+1;j<=num_ctg;j++)
			if(length_array[j]!=length_array[i])
				break;
		length_array[++diff_len] = length_array[i];
		flag_array[diff_len] = i;
		i = j - 1;
	}
	/*
	for(i=1;i<=num_ctg;i++)
		flag_array[i] = 0;
	*/
	contig_array = (CONTIG *)ckalloc((num_ctg+1)*sizeof(CONTIG));

	//load edges
	index = 0;
	rewind(fp);
	while(fgets(line,sizeof(line),fp)!=NULL){
		if(line[0]=='>'){
//			if(overlaplen<=31)
//				sscanf(line,">length %u,%llx,%llx,%d,%d",&length,&from_kmer,&to_kmer,&bal_ed,&cvg);
//			else
				sscanf(line,">length %u,%d,%d",&length,&bal_ed,&cvg);
			newIndex = uniqueLenSearch(length_array,flag_array,diff_len,length);
			index_array[++index]=newIndex;
			
			contig_array[newIndex].length = length;
			contig_array[newIndex].bal_edge = bal_ed + 1;
			contig_array[newIndex].downwardConnect = NULL;
			contig_array[newIndex].mask = 0;
			contig_array[newIndex].flag = 0;
			contig_array[newIndex].arcs = NULL;
			contig_array[newIndex].seq = NULL;
			contig_array[newIndex].multi = 0;
			contig_array[newIndex].inSubGraph = 0;
			contig_array[newIndex].cvg = cvg/10; 
			if(cvg){
				counter += length;
				cvgSum += cvg*length;
			}
			fprintf(out_fp,"%d %d %d\n",index,newIndex,contig_array[newIndex].bal_edge);
		}
	}
	if(counter)
		//cvgAvg = cvgSum/counter > 2 ? cvgSum/counter : 3;
		cvgAvg = cvgSum/counter/10 > 2 ? cvgSum/counter/10 : 3;
	
	//mark repeats
	int bal_i;
	/*if(maskRep){
		counter = 0;
		for(i=1;i<=num_ctg;i++){
			bal_i = getTwinCtg(i);
			if((contig_array[i].cvg+contig_array[bal_i].cvg)>4*cvgAvg){
				contig_array[i].mask = 1;
				contig_array[bal_i].mask = 1;
				counter += 2;
			}
			if(isSmallerThanTwin(i))
				i++;
		}
		printf("average contig coverage : %d. Number of contig(s) masked because of high coverage: %llx\n",
						cvgAvg,counter);
	}*/

	counter = 0;
	for(i=1;i<=num_ctg;i++){
		if(contig_array[i].mask)
			continue;
		bal_i = getTwinCtg(i);
		if(contig_array[i].length<ctg_short){
			contig_array[i].mask = 1;
			contig_array[bal_i].mask = 1;
			counter += 2;
		}
		if(isSmallerThanTwin(i))
			i++;
	}
	//printf("Mask contigs shorter than %d, %lld contig masked\n",
		//					ctg_short,counter);
		printf("[%s]Number of contig(s) masked because of shortie: %lld\n",__FUNCTION__,counter);
	//loadArcs(graphfile);
	//tipsCount();
	loadContig(graphfile);
	//printf("done loading updated edges\n");
	fflush(stdout);
	free((void *)length_array);
	free((void *)flag_array);
	fclose(fp);
	fclose(out_fp);
}

/*static void add1Arc(unsigned int from_ed,unsigned int to_ed, unsigned int weight)
{
	preARC *parc;
	unsigned int from_c = index_array[from_ed];
	unsigned int to_c = index_array[to_ed];
	
	parc = allocatePreArc(to_c);
	parc->multiplicity = weight;
	parc->next = contig_array[from_c].arcs;
	contig_array[from_c].arcs = parc;
}*/

/*void loadArcs(char *graphfile)
{
	FILE *fp;
	char name[256],line[1024];
	unsigned int target,weight;
	unsigned int from_ed;
	char *seg;

	sprintf(name,"%s.Arc",graphfile);
	fp = ckopen(name,"r");

	createPreArcMemManager();
	arcCounter = 0;
	while(fgets(line,sizeof(line),fp)!=NULL){
		seg = strtok(line," ");
		from_ed = atoi(seg);
		//printf("%d\n",from_ed);
		while((seg=strtok(NULL," "))!=NULL){
			target = atoi(seg);
			seg = strtok(NULL," ");
			weight = atoi(seg);
			add1Arc(from_ed,target,weight);
			
		}
	}
	printf("%lld arcs loaded\n",arcCounter);
	fclose(fp);
}*/

void loadContig(char *graphfile)
{
	//fprintf(stderr,"[%s]entering this function\n",__FUNCTION__);
	char c,name[256],line[1024],*tightSeq=NULL;
	FILE *fp;
	int n=0,length,index=-1,edgeno;
	unsigned int i;
	unsigned int newIndex;
	
	sprintf(name,"%s.contig",graphfile);
	fp = ckopen(name,"r");

	while(fgets(line,sizeof(line),fp)!=NULL){
		if(line[0]=='>'){
			if(index>=0){
				newIndex = index_array[edgeno];
				contig_array[newIndex].seq = tightSeq;
			}
			n=0;
			index++;
			sscanf(line+1,"%d %s %d",&edgeno,name,&length);
			//printf("contig %d, length %d\n",edgeno,length);
			tightSeq = (char *)ckalloc((length/4+1)*sizeof(char));
			//fprintf(stderr,"[%s]loaded %d.\n",__FUNCTION__,edgeno);
		}else{
			int tmp_len=strlen(line);
			for(i=0;i<tmp_len;i++){
				if(line[i]>='a' && line[i]<='z'){
					c = base2int(line[i]-'a'+'A');
					writeChar2tightString(c,tightSeq,n++);
				}
				else if(line[i]>='A' && line[i]<='Z'){
					c = base2int(line[i]);
					writeChar2tightString(c,tightSeq,n++);
				}
			}
		}

	}
	if(index>=0){
		newIndex = index_array[edgeno];
		contig_array[newIndex].seq = tightSeq;
	}
	printf("[%s]input %d contigs\n",__FUNCTION__,index+1);
	fclose(fp);
	
	//printf("the %dth contig with index 107\n",index);
}
void freeContig_array()
{
	if(!contig_array)
		return;

	unsigned int i;

	for(i=1;i<=num_ctg;i++){
		if(contig_array[i].seq)
			free((void *)contig_array[i].seq);
		if(contig_array[i].closeReads)
			freeStack(contig_array[i].closeReads);
	}

	free((void *)contig_array);
	contig_array = NULL;
}
/*
void loadCvg(char *graphfile)
{
	char name[256],line[1024];
	FILE *fp;
	int cvg;
	unsigned int newIndex,edgeno,bal_ctg;
	
	sprintf(name,"%s.contigCVG",graphfile);
	fp = fopen(name,"r");
	if(!fp){
		printf("contig coverage file %s is not found!\n",name);
		return;
	}

	while(fgets(line,sizeof(line),fp)!=NULL){
		if(line[0]=='>'){
			sscanf(line+1,"%d %d",&edgeno,&cvg);
			newIndex = index_array[edgeno];
			cvg = cvg <= 255 ? cvg:255;
			contig_array[newIndex].multi = cvg;
			bal_ctg = getTwinCtg(newIndex);
			contig_array[bal_ctg].multi= cvg;
		}
	}
	fclose(fp);
}
*/

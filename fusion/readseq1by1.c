#include "stdinc.h"
#include "newhash.h"
#include "extfunc.h"
#include "extvab.h"

static char src_rc_seq[1024];
extern long long single_count;
extern long long single_map;
void readseq1by1(char *src_seq, char *src_name, int *len_seq, FILE *fp,long long num_seq)
{
	int	i,k, n,strL;
	char	c;
	char	str[5000];

	n = 0;
	k = num_seq;
	while(fgets(str, 4950, fp))	{
		if(str[0] == '#')	continue;
		if(str[0] == '>')	{
			/*
			if(k >= 0)	{  // if this isn't the first '>' in the file
				*len_seq = n;
			}
			*/
			*len_seq = n;
			n = 0;
			sscanf(&str[1],"%s",src_name);
			return;
		} else {
			strL = strlen(str);
			if(strL+n>maxReadLen)
				strL = maxReadLen - n;
			for(i = 0; i < strL; i ++)	{
				if(str[i] >= 'a' && str[i] <= 'z') {
					c = base2int(str[i]-'a'+'A');
					src_seq[n ++] = c;
				} else if(str[i] >= 'A' && str[i] <= 'Z') {
					c = base2int(str[i]);
					src_seq[n ++] = c;
				          // after pre-process all the symbles would be a,g,c,t,n in lower or upper case.
				} else if(str[i]=='.') {
					c = base2int('A');
					src_seq[n ++] = c;
				}          // after pre-process all the symbles would be a,g,c,t,n in lower or upper case.
			}
			//printf("%d: %d\n",k,n);
		}
	}

	if(k >= 0){
		*len_seq = n;
		return;
	}
	*len_seq = 0;
}


void read_one_sequence(FILE *fp, long long *T, char **X)

{

	char *fasta,*src_name; //point to fasta array
	int num_seq,len,name_len,min_len;

	num_seq = readseqpar(&len,&min_len,&name_len,fp);
	if(num_seq<1){
		printf("no fasta sequence in file\n");
		*T = 0;
		return;
	}
	fasta = (char *)ckalloc(len*sizeof(char));
	src_name = (char *)ckalloc((name_len+1)*sizeof(char));
	rewind(fp);

	readseq1by1(fasta,src_name,&len,fp,-1);
	readseq1by1(fasta,src_name,&len,fp,0);
	
	*X = fasta;
	*T = len;
	free((void *)src_name);
}

long long multiFileParse(int *max_leg, int *min_leg,int *max_name_leg, FILE *fp)
{
	
	char str[5000];
	FILE *freads;
	int slen;
	long long counter = 0;
	*max_name_leg = *max_leg = 1;
	*min_leg = 1000;
	while(fgets(str,4950,fp)){
		slen = strlen(str);
		str[slen-1] = str[slen];
		freads = ckopen(str,"r");
		counter += readseqpar(max_leg,min_leg,max_name_leg,freads);
		fclose(freads);
	}
	return counter;
}

long long readseqpar(int *max_leg, int *min_leg,int *max_name_leg, FILE *fp)
{
	int	l, n;
	long long k;
	char	str[5000], src_name[5000];


	n = 0;
	k = -1;
	while(fgets(str, 4950, fp))	{
		if(str[0] == '>')	{
			if(k >= 0)	{
				if(n > *max_leg)	
					*max_leg = n;
				if(n < *min_leg)	
					*min_leg = n;
			
			}
			n = 0;
			k ++;
			sscanf(&str[1], "%s", src_name);
			if((l = strlen(src_name)) > *max_name_leg)		
				*max_name_leg = l;
		} else {
			n += strlen(str)-1;
		}
	}
	if(n > *max_leg)	
		*max_leg = n;
	
	if(n < *min_leg)	
		*min_leg = n;

	k ++;
	return(k);
}

void read1seqfq(char *src_seq, char *src_name, int *len_seq, FILE *fp)
{
	int	i,n,strL;
	char	c;
	char	str[5000];
	boolean flag=0;

	while(fgets(str, 4950, fp))	{
		if(str[0]=='@'){
			flag = 1;
			sscanf(&str[1],"%s",src_name);
			break;
		}
	}

	if(!flag){   //last time reading fq file get this
		*len_seq = 0;
		return;
	}

	n = 0;
	while(fgets(str, 4950, fp)){
		if(str[0] == '+')	{
			fgets(str,4950,fp);   // pass quality value line
			*len_seq = n;
			return;
		} else {
			strL = strlen(str);
			if(strL+n>maxReadLen)
				strL = maxReadLen - n;
			for(i = 0; i < strL; i ++)	{
				if(str[i] >= 'a' && str[i] <= 'z') {
					c = base2int(str[i]-'a'+'A');
					src_seq[n ++] = c;
				} else if(str[i] >= 'A' && str[i] <= 'Z') {
					c = base2int(str[i]);
					src_seq[n ++] = c;
				          // after pre-process all the symbles would be a,g,c,t,n in lower or upper case.
				} else if(str[i]=='.') {
					c = base2int('A');
					src_seq[n ++] = c;
				}          // after pre-process all the symbles would be a,g,c,t,n in lower or upper case.
			}
			//printf("%d: %d\n",k,n);
		}
	}

	*len_seq = n;
	return;
}

// find the next file to open in libs
static int nextValidIndex(int libNo,boolean pair,unsigned char asm_ctg)
{
	int i=libNo;
	
	while(i<num_libs){
		if(asm_ctg==1&&(lib_array[i].asm_flag!=1&&lib_array[i].asm_flag!=3)){
			i++;
			continue;
		}else if(asm_ctg==0&&(lib_array[i].asm_flag!=2&&lib_array[i].asm_flag!=3)){
			i++;
			continue;
		}else if(asm_ctg>1&&lib_array[i].asm_flag!=asm_ctg){  // reads for other purpose 
			i++;
			continue;
		}
		if(lib_array[i].curr_type==1&&
			lib_array[i].curr_index<lib_array[i].num_a1_file)
			return i;
		if(lib_array[i].curr_type==2&&
			lib_array[i].curr_index<lib_array[i].num_q1_file)
			return i;
		if(lib_array[i].curr_type==3&&
			lib_array[i].curr_index<lib_array[i].num_p_file)
			return i;
		if(pair){
			if(lib_array[i].curr_type<3){
				lib_array[i].curr_type++;
				lib_array[i].curr_index = 0;
			}else
				i++;
			continue;
		}
		if(lib_array[i].curr_type==4&&
			lib_array[i].curr_index<lib_array[i].num_s_a_file)
			return i;
		if(lib_array[i].curr_type==5&&
			lib_array[i].curr_index<lib_array[i].num_s_q_file)
			return i;
		if(lib_array[i].curr_type<5){
			lib_array[i].curr_type++;
			lib_array[i].curr_index = 0;
		}else
			i++;
	}//for each lib
		
	return i;
}

static FILE *openFile4read(char *fname)
{
	FILE *fp;
	if(strlen(fname)>3&&strcmp(fname+strlen(fname)-3,".gz")==0){
		char *cmd = (char *)ckalloc((strlen(fname)+20)*sizeof(char));
		sprintf(cmd,"gzip -dc %s",fname);
		fp = popen(cmd,"r");
		free(cmd);
		return fp;
	}else{
		return ckopen(fname,"r");	
	}
	
}

void openFileInLib(int libNo)
{
	int i = libNo;
	if(lib_array[i].curr_type==1){
		printf("[%s]opened file:\n %s\n",
					__FUNCTION__,lib_array[i].a1_fname[lib_array[i].curr_index]);
		printf("[%s]opened file:\n %s\n",
					__FUNCTION__,lib_array[i].a2_fname[lib_array[i].curr_index]);
		lib_array[i].fp1 = openFile4read(lib_array[i].a1_fname[lib_array[i].curr_index]);	
		lib_array[i].fp2 = openFile4read(lib_array[i].a2_fname[lib_array[i].curr_index]);	
		lib_array[i].curr_index++;
		lib_array[i].paired = 1;
	}else if(lib_array[i].curr_type==2){
		printf("[%s]opened file:\n %s\n",
					__FUNCTION__,lib_array[i].q1_fname[lib_array[i].curr_index]);
		printf("[%s]opened file:\n %s\n",
					__FUNCTION__,lib_array[i].q2_fname[lib_array[i].curr_index]);
		lib_array[i].fp1 = openFile4read(lib_array[i].q1_fname[lib_array[i].curr_index]);	
		lib_array[i].fp2 = openFile4read(lib_array[i].q2_fname[lib_array[i].curr_index]);	
		lib_array[i].curr_index++;
		lib_array[i].paired = 1;
	}else if(lib_array[i].curr_type==3){
		printf("[%s]opened file:\n %s\n",
					lib_array[i].p_fname[lib_array[i].curr_index]);
		lib_array[i].fp1 = openFile4read(lib_array[i].p_fname[lib_array[i].curr_index]);	
		lib_array[i].curr_index++;
		lib_array[i].paired = 0;
	}else if(lib_array[i].curr_type==4){
		printf("[%s]opened file:\n %s\n",
					__FUNCTION__,lib_array[i].s_a_fname[lib_array[i].curr_index]);
		lib_array[i].fp1 = openFile4read(lib_array[i].s_a_fname[lib_array[i].curr_index]);	
		lib_array[i].curr_index++;
		lib_array[i].paired = 0;
	}else if(lib_array[i].curr_type==5){
		printf("[%s]opened file:\n %s\n",
					__FUNCTION__,lib_array[i].s_q_fname[lib_array[i].curr_index]);
		lib_array[i].fp1 = openFile4read(lib_array[i].s_q_fname[lib_array[i].curr_index]);	
		lib_array[i].curr_index++;
		lib_array[i].paired = 0;
	}
	
}

static void reverse2k(char *src_seq,int len_seq)
{
	if(!len_seq)
		return;
	
	int i;
	reverseComplementSeq(src_seq,len_seq,src_rc_seq);

	for(i=0;i<len_seq;i++)
		src_seq[i] = src_rc_seq[i];
}

static void closeFp1InLab(int libNo)
{
	int ftype = lib_array[libNo].curr_type;
	int index = lib_array[libNo].curr_index-1;
	char *fname;
	if(ftype==1)
		fname = lib_array[libNo].a1_fname[index];
	else if(ftype==2)
		fname = lib_array[libNo].q1_fname[index];
	else if(ftype==3)
		fname = lib_array[libNo].p_fname[index];
	else if(ftype==4)
		fname = lib_array[libNo].s_a_fname[index];
	else if(ftype==5)
		fname = lib_array[libNo].s_q_fname[index];
	else
		return;
	if(strlen(fname)>3&&strcmp(fname+strlen(fname)-3,".gz")==0)
		pclose(lib_array[libNo].fp1);
	else
		fclose(lib_array[libNo].fp1);
}

static void closeFp2InLab(int libNo)
{
	int ftype = lib_array[libNo].curr_type;
	int index = lib_array[libNo].curr_index-1;
	char *fname;
	if(ftype==1)
		fname = lib_array[libNo].a2_fname[index];
	else if(ftype==2)
		fname = lib_array[libNo].q2_fname[index];
	else
		return;
	if(strlen(fname)>3&&strcmp(fname+strlen(fname)-3,".gz")==0)
		pclose(lib_array[libNo].fp2);
	else
		fclose(lib_array[libNo].fp2);
}

boolean read1seqInLib(char *src_seq, char *src_name, int *len_seq, int *libNo,boolean pair,unsigned char asm_ctg)
{
	int i = *libNo;
	int prevLib = i;

	if(!lib_array[i].fp1  // file1 does not exist 
		||(lib_array[i].curr_type!=1&&feof(lib_array[i].fp1))   // file1 reaches end and not type1
			||(lib_array[i].curr_type==1&&feof(lib_array[i].fp1)&&feof(lib_array[i].fp2))){//f1&f2 reaches end
		if(lib_array[i].fp1&&feof(lib_array[i].fp1)){
				 closeFp1InLab(i);
				 //printf("[%s]%d reads in current file , (%.1f) map-rate .\n",__FUNCTION__,single_count,single_map/single_count);
				 single_count=single_map=0;
		}
		if(lib_array[i].fp2&&feof(lib_array[i].fp2)){
				closeFp2InLab(i);
				//printf("[%s]%d reads in current file , (%.1f) map-rate .\n",__FUNCTION__,single_count,single_map/single_count);
				single_count=single_map=0;
		}
		
		*libNo = nextValidIndex(i,pair,asm_ctg);
		i = *libNo;
		if(lib_array[i].rd_len_cutoff>0)
			maxReadLen = lib_array[i].rd_len_cutoff<maxReadLen4all ?
								lib_array[i].rd_len_cutoff:maxReadLen4all;
		else
			maxReadLen = maxReadLen4all;

		//record insert size info
		//printf("from lib %d to %d, read %lld to %ld\n",prevLib,i,readNumBack,n_solexa);
		if(pair&&i!=prevLib){
			if(readNumBack<n_solexa){
				pes[gradsCounter].PE_bound = n_solexa;
				pes[gradsCounter].rank = lib_array[prevLib].rank;
				pes[gradsCounter].pair_num_cut = lib_array[prevLib].pair_num_cut;
				pes[gradsCounter++].insertS = lib_array[prevLib].avg_ins;
				readNumBack = n_solexa;
			}		
		}
		if(i>=num_libs)
			return 0;
		openFileInLib(i);

		if(lib_array[i].curr_type==1){
			readseq1by1(src_seq, src_name,len_seq, lib_array[i].fp1,-1);
			readseq1by1(src_seq, src_name,len_seq, lib_array[i].fp2,-1);
		}else if(lib_array[i].curr_type==3||lib_array[i].curr_type==4)
			readseq1by1(src_seq, src_name,len_seq, lib_array[i].fp1,-1);

	}
	if(lib_array[i].curr_type==1){
		if(lib_array[i].paired==1){
			readseq1by1(src_seq, src_name,len_seq, lib_array[i].fp1,1);
			if(lib_array[i].reverse)
				reverse2k(src_seq,*len_seq);
			lib_array[i].paired = 2;
			if(*len_seq>0||!feof(lib_array[i].fp1)){
				n_solexa++;
				return 1;
			}
			else
				return read1seqInLib(src_seq,src_name,len_seq,libNo,pair,asm_ctg);
		}else{
			readseq1by1(src_seq, src_name,len_seq, lib_array[i].fp2,1);

			if(lib_array[i].reverse)
				reverse2k(src_seq,*len_seq);
			lib_array[i].paired = 1;
			n_solexa++;
			return 1; //can't fail to read a read2
		}
	}
	if(lib_array[i].curr_type==2){
		if(lib_array[i].paired==1){
			read1seqfq(src_seq, src_name,len_seq, lib_array[i].fp1);
			/*
			if(*len_seq>0){
				for(j=0;j<*len_seq;j++)
					printf("%c",int2base(src_seq[j]));	
				printf("\n");
			}
			*/
			if(lib_array[i].reverse)
				reverse2k(src_seq,*len_seq);
			lib_array[i].paired = 2;
			if(*len_seq>0||!feof(lib_array[i].fp1)){
				n_solexa++;
				return 1;
			}else
				return read1seqInLib(src_seq,src_name,len_seq,libNo,pair,asm_ctg);
		}else{
			read1seqfq(src_seq, src_name,len_seq, lib_array[i].fp2);
			if(lib_array[i].reverse)
				reverse2k(src_seq,*len_seq);
			lib_array[i].paired = 1;
			n_solexa++;
			return 1; //can't fail to read a read2
		}
	}
	if(lib_array[i].curr_type==5)
		read1seqfq(src_seq, src_name,len_seq, lib_array[i].fp1);
	else{
		readseq1by1(src_seq, src_name,len_seq, lib_array[i].fp1,1);
	}
	/*
	int t;
	for(t=0;t<*len_seq;t++)
		printf("%d",src_seq[t]);
	printf("\n");
	*/
	if(lib_array[i].reverse)
		reverse2k(src_seq,*len_seq);
	if(*len_seq>0||!feof(lib_array[i].fp1)){
		n_solexa++;
		return 1;
	}else
		return read1seqInLib(src_seq,src_name,len_seq,libNo,pair,asm_ctg);
}

#include "stdinc.h"
#include "newhash.h"
#include "extfunc.h"
#include "extvab.h"

static char tabs[2][1024];

static boolean splitColumn(char *line)
{
	int len = strlen(line);
	int i=0,j;
	int tabs_n = 0;

	while(i<len){
		if(line[i]>=32&&line[i]<=126&&line[i]!='='){
			j=0;
			while(i<len&&line[i]>=32&&line[i]<=126&&line[i]!='='){
				tabs[tabs_n][j++] = line[i];	
				i++;
			}
			tabs[tabs_n][j] = '\0';	
			tabs_n++;
			if(tabs_n==2)
				return 1;
		}
		i++;
	}
	if(tabs_n==2)
		return 1;
	else
		return 0;
}

static int cmp_lib(const void *a,const void *b)
{
	LIB_INFO *A,*B;
	A = (LIB_INFO *)a;
	B = (LIB_INFO *)b;

	if(A->avg_ins>B->avg_ins)
		return 1;
	else if(A->avg_ins==B->avg_ins)
		return 0;
	else
		return -1;
}

void scan_libInfo(char *libfile)
{
	FILE *fp;
	char line[1024],ch;
	int i,j,index;
	int libCounter;
	boolean flag;
	
	fp = ckopen(libfile,"r");
	num_libs = 0;
	while(fgets(line,1024,fp)){
		ch = line[5];
		line[5] = '\0';
		if(strcmp(line,"[LIB]")==0)
			num_libs++;
		if(!num_libs){
			line[5] = ch;
			flag = splitColumn(line);
			if(!flag)
				continue;
			if(strcmp(tabs[0],"max_rd_len")==0)
				maxReadLen = atoi(tabs[1]);
		}
	}
//count file numbers of each type
	lib_array = (LIB_INFO *)ckalloc(num_libs*sizeof(LIB_INFO));
	for(i=0;i<num_libs;i++){
		lib_array[i].asm_flag = 3;
		lib_array[i].rd_len_cutoff = 0;
		lib_array[i].rank = 0;
		lib_array[i].pair_num_cut = 1;
		lib_array[i].map_len = 0;
		lib_array[i].num_s_a_file=0;	
		lib_array[i].num_s_q_file=0;	
		lib_array[i].num_p_file=0;	
		
		lib_array[i].num_a1_file=0;	
		lib_array[i].num_a2_file=0;	
		lib_array[i].num_q1_file=0;	
		lib_array[i].num_q2_file=0;	
	}
	libCounter = -1;
	rewind(fp);
	
	i = -1;
	while(fgets(line,1024,fp)){
		ch = line[5];
		line[5] = '\0';
		if(strcmp(line,"[LIB]")==0){
			i++;
			continue;
		}
		line[5] = ch;
		flag = splitColumn(line);
		if(!flag)
			continue;

		if(strcmp(tabs[0],"f1")==0){
			lib_array[i].num_a1_file++;
		}else if(strcmp(tabs[0],"q1")==0){
			lib_array[i].num_q1_file++;
		}else if(strcmp(tabs[0],"f2")==0){
			lib_array[i].num_a2_file++;
		}else if(strcmp(tabs[0],"q2")==0){
			lib_array[i].num_q2_file++;
		}else if(strcmp(tabs[0],"f")==0){
			lib_array[i].num_s_a_file++;
		}else if(strcmp(tabs[0],"q")==0){
			lib_array[i].num_s_q_file++;
		}else if(strcmp(tabs[0],"p")==0){
			lib_array[i].num_p_file++;
		}
	}
//allocate memory for filenames
	for(i=0;i<num_libs;i++){
		
		if(lib_array[i].num_s_a_file){
			lib_array[i].s_a_fname = (char **)ckalloc(lib_array[i].num_s_a_file*sizeof(char *));
			for(j=0;j<lib_array[i].num_s_a_file;j++)
				lib_array[i].s_a_fname[j] = (char *)ckalloc(1024*sizeof(char));
		}
		if(lib_array[i].num_s_q_file){
			lib_array[i].s_q_fname = (char **)ckalloc(lib_array[i].num_s_q_file*sizeof(char *));
			for(j=0;j<lib_array[i].num_s_q_file;j++)
				lib_array[i].s_q_fname[j] = (char *)ckalloc(1024*sizeof(char));
		}
		if(lib_array[i].num_p_file){
			lib_array[i].p_fname = (char **)ckalloc(lib_array[i].num_p_file*sizeof(char *));
			for(j=0;j<lib_array[i].num_p_file;j++)
				lib_array[i].p_fname[j] = (char *)ckalloc(1024*sizeof(char));
		}
		if(lib_array[i].num_a1_file){
			lib_array[i].a1_fname = (char **)ckalloc(lib_array[i].num_a1_file*sizeof(char *));
			for(j=0;j<lib_array[i].num_a1_file;j++)
				lib_array[i].a1_fname[j] = (char *)ckalloc(1024*sizeof(char));
		}
		if(lib_array[i].num_a2_file){
			lib_array[i].a2_fname = (char **)ckalloc(lib_array[i].num_a2_file*sizeof(char *));
			for(j=0;j<lib_array[i].num_a2_file;j++)
				lib_array[i].a2_fname[j] = (char *)ckalloc(1024*sizeof(char));
		}
		if(lib_array[i].num_q1_file){
			lib_array[i].q1_fname = (char **)ckalloc(lib_array[i].num_q1_file*sizeof(char *));
			for(j=0;j<lib_array[i].num_q1_file;j++)
				lib_array[i].q1_fname[j] = (char *)ckalloc(1024*sizeof(char));
		}
		if(lib_array[i].num_q2_file){
			lib_array[i].q2_fname = (char **)ckalloc(lib_array[i].num_q2_file*sizeof(char *));
			for(j=0;j<lib_array[i].num_q2_file;j++)
				lib_array[i].q2_fname[j] = (char *)ckalloc(1024*sizeof(char));
		}
	}

// get file names
	for(i=0;i<num_libs;i++){
		lib_array[i].curr_type=1;
		lib_array[i].curr_index=0;
		lib_array[i].fp1=NULL;
		lib_array[i].fp2=NULL;

		lib_array[i].num_s_a_file=0;	
		lib_array[i].num_s_q_file=0;	
		lib_array[i].num_p_file=0;	
		
		lib_array[i].num_a1_file=0;	
		lib_array[i].num_a2_file=0;	
		lib_array[i].num_q1_file=0;	
		lib_array[i].num_q2_file=0;	
	}
	libCounter = -1;
	rewind(fp);
	
	i = -1;
	while(fgets(line,1024,fp)){
		ch = line[5];
		line[5] = '\0';
		if(strcmp(line,"[LIB]")==0){
			i++;
			continue;
		}
		line[5] = ch;
		flag = splitColumn(line);
		if(!flag)
			continue;

		if(strcmp(tabs[0],"f1")==0){
			index = lib_array[i].num_a1_file++;
			strcpy(lib_array[i].a1_fname[index],tabs[1]);
		}else if(strcmp(tabs[0],"q1")==0){
			index = lib_array[i].num_q1_file++;
			strcpy(lib_array[i].q1_fname[index],tabs[1]);
		}else if(strcmp(tabs[0],"f2")==0){
			index = lib_array[i].num_a2_file++;
			strcpy(lib_array[i].a2_fname[index],tabs[1]);
		}else if(strcmp(tabs[0],"q2")==0){
			index = lib_array[i].num_q2_file++;
			strcpy(lib_array[i].q2_fname[index],tabs[1]);
		}else if(strcmp(tabs[0],"f")==0){
			index = lib_array[i].num_s_a_file++;
			strcpy(lib_array[i].s_a_fname[index],tabs[1]);
		}else if(strcmp(tabs[0],"q")==0){
			index = lib_array[i].num_s_q_file++;
			strcpy(lib_array[i].s_q_fname[index],tabs[1]);
		}else if(strcmp(tabs[0],"p")==0){
			index = lib_array[i].num_p_file++;
			strcpy(lib_array[i].p_fname[index],tabs[1]);
		}else if(strcmp(tabs[0],"min_ins")==0)
			lib_array[i].min_ins = atoi(tabs[1]);
		else if(strcmp(tabs[0],"max_ins")==0)
			lib_array[i].max_ins = atoi(tabs[1]);
		else if(strcmp(tabs[0],"avg_ins")==0)
			lib_array[i].avg_ins = atoi(tabs[1]);
		else if(strcmp(tabs[0],"rd_len_cutoff")==0)
			lib_array[i].rd_len_cutoff = atoi(tabs[1]);
		else if(strcmp(tabs[0],"reverse_seq")==0)
			lib_array[i].reverse = atoi(tabs[1]);
		else if(strcmp(tabs[0],"asm_flags")==0)
			lib_array[i].asm_flag = atoi(tabs[1]);
		else if(strcmp(tabs[0],"rank")==0)
			lib_array[i].rank = atoi(tabs[1]);
		else if(strcmp(tabs[0],"pair_num_cutoff")==0)
			lib_array[i].pair_num_cut = atoi(tabs[1]);
		else if(strcmp(tabs[0],"rd_len_cutoff")==0)
			lib_array[i].rd_len_cutoff = atoi(tabs[1]);
		else if(strcmp(tabs[0],"map_len")==0)
			lib_array[i].map_len = atoi(tabs[1]);
		
	}

	fclose(fp);

	qsort(&lib_array[0],num_libs,sizeof(LIB_INFO),cmp_lib);
}

int getMaxLongReadLen(int num_libs)
{
	int i;
	int maxLong=0;
	boolean Has=0;
	for(i=0;i<num_libs;i++){
		if(lib_array[i].asm_flag!=4)
			continue;
		Has = 1;
		maxLong = maxLong<lib_array[i].rd_len_cutoff ? lib_array[i].rd_len_cutoff:maxLong;
	}
	if(!Has)
		return maxLong;
	else
		return maxLong>0 ? maxLong:maxReadLen;
}

void free_libs()
{

	if(!lib_array)
		return;

	int i,j;
	for(i=0;i<num_libs;i++){
		//printf("[LIB] %d, avg_ins %d, reverse %d \n",i,lib_array[i].avg_ins,lib_array[i].reverse);
		if(lib_array[i].num_s_a_file){
			//printf("%d single fasta files\n",lib_array[i].num_s_a_file);
			for(j=0;j<lib_array[i].num_s_a_file;j++)
				free((void *)lib_array[i].s_a_fname[j]);
			free((void *)lib_array[i].s_a_fname);
		}
		if(lib_array[i].num_s_q_file){
			//printf("%d single fastq files\n",lib_array[i].num_s_q_file);
			for(j=0;j<lib_array[i].num_s_q_file;j++)
				free((void *)lib_array[i].s_q_fname[j]);
			free((void *)lib_array[i].s_q_fname);
		}
		if(lib_array[i].num_p_file){
			//printf("%d paired fasta files\n",lib_array[i].num_p_file);
			for(j=0;j<lib_array[i].num_p_file;j++)
				free((void *)lib_array[i].p_fname[j]);
			free((void *)lib_array[i].p_fname);
		}
		if(lib_array[i].num_a1_file){
			//printf("%d read1 fasta files\n",lib_array[i].num_a1_file);
			for(j=0;j<lib_array[i].num_a1_file;j++)
				free((void *)lib_array[i].a1_fname[j]);
			free((void *)lib_array[i].a1_fname);
		}
		if(lib_array[i].num_a2_file){
			//printf("%d read2 fasta files\n",lib_array[i].num_a2_file);
			for(j=0;j<lib_array[i].num_a2_file;j++)
				free((void *)lib_array[i].a2_fname[j]);
			free((void *)lib_array[i].a2_fname);
		}
		if(lib_array[i].num_q1_file){
			//printf("%d read1 fastq files\n",lib_array[i].num_q1_file);
			for(j=0;j<lib_array[i].num_q1_file;j++)
				free((void *)lib_array[i].q1_fname[j]);
			free((void *)lib_array[i].q1_fname);
		}
		if(lib_array[i].num_q2_file){
			//printf("%d read2 fastq files\n",lib_array[i].num_q2_file);
			for(j=0;j<lib_array[i].num_q2_file;j++)
				free((void *)lib_array[i].q2_fname[j]);
			free((void *)lib_array[i].q2_fname);
		}
	}
	
	num_libs = 0;
	free((void *)lib_array);
}

void alloc_pe_mem(int gradsCounter)
{
	if(gradsCounter)
		pes = (PE_INFO *)ckalloc(gradsCounter*sizeof(PE_INFO));
}

void free_pe_mem()
{
	if(pes){
		free((void *)pes);
		pes = NULL;
	}
}


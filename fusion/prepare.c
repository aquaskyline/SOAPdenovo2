#include "stdinc.h" 
#include "newhash.h"
#include "extfunc.h" 
#include "extvab.h" 
#include "ctype.h"
boolean upper_rev(char *in,int in_len);
void print_seq(FILE *out_file,char *sequence , int sequence_len);
char rev[]={0,0,0,0,0,0,0,0,0,0,//0
					0,0,0,0,0,0,0,0,0,0,//10
					0,0,0,0,0,0,0,0,0,0,//20
					0,0,0,0,0,0,0,0,0,0,//30
					0,0,0,0,0,0,0,0,0,0,//40
					0,0,0,0,0,0,0,0,0,0,//50
					0,0,0,0,0,'T',0,'G',0,0,//60
					0,'C',0,0,0,0,0,0,'N',0,//70
					0,0,0,0,'A',0,0,0,0,0,//80
					0,0,0,0,0,0,0,0,0,0,//90
					0,0,0,0,0,0,0,0,0,0,};//100
typedef struct io_ctg{
	char *seq;
	int len;
	int bal;
	char *name;
}IO_CTG;

static int cmp_ctg(const void *a,const void *b){
	IO_CTG *A=(IO_CTG *)a;
	IO_CTG *B=(IO_CTG *)b;
	return A->len-B->len;
}

int data_prepare(){
	char file_name[256];

	FILE *basic;
	sprintf(file_name,"%s.preGraphBasic",graphfile);
	basic=ckopen(file_name,"w");
	fprintf(basic,"VERTEX 605681 K %d",overlaplen);
	fprintf(basic,"\nEDGEs 1861091\n\nMaxReadLen 100 MinReadLen 0 MaxNameLen 256\n");
	fclose(basic);
	
	//char **ctg_seq=(char **)ckalloc(100000000*sizeof(char *));
	//int *ctg_bal=(int *)ckalloc(100000000*sizeof(int));
	//int *ctg_len=(int *)ckalloc(100000000*sizeof(int));

	FILE *ctg_fp;
	ctg_fp=ckopen(ctg_file,"r");
	FILE *update,*index,*new_ctg;
	sprintf(file_name, "%s.contig", graphfile);
	new_ctg=ckopen(file_name,"w");
	FILE *conver;
	sprintf(file_name,"%s.conver", graphfile);
	conver=ckopen(file_name,"w");
	
	
	char *line;
	line= (char *)ckalloc(100000000*sizeof(char ));
	char orig_name[1024];
	char *seq;
	IO_CTG *pre_ctg=(IO_CTG *)ckalloc(1000000000*sizeof(IO_CTG));

	seq=(char *)malloc(1000000000*sizeof(char));
	int cul_id=1;
	int total=0;
	fgets(line,100000000*sizeof(char ),ctg_fp);
	sscanf(line,">%s",orig_name);
	int len=0;
	//fprintf(stderr,"reach here %d\n",__LINE__);
	while(fgets(line,100000000*sizeof(char ),ctg_fp)!=NULL){
		if(line[0]=='>'){
			if(len<overlaplen){
				sscanf(line,">%s",orig_name);
				seq[0]='\0';
				len=0;
				continue;
			}
			
			boolean flag=upper_rev(seq,len);
			//fprintf(new_ctg,">%d length %d\n",cul_id,len);
			//fprintf(conver,"%s\t%d\t%d\n",orig_name,cul_id,len);
			//print_seq(new_ctg,seq,len);
			//fprintf(new_ctg,"%s\n",seq);
			char *one_seq=(char *)ckalloc((len+100)*sizeof(char));
			strcpy(one_seq,seq);
			
			if(flag==0){
				pre_ctg[++total].seq=one_seq;
				pre_ctg[total].bal=2;
				pre_ctg[total].len=len;
				pre_ctg[total].name=(char *)malloc((strlen(orig_name)+1)*sizeof(char));
				strcpy(pre_ctg[total].name,orig_name);
				//pre_ctg[++cul_id].bal=0;
				cul_id+=2;		
			}else{
				pre_ctg[++total].seq=one_seq;
				pre_ctg[total].len=len;
				pre_ctg[total].bal=1;
				pre_ctg[total].name=(char *)malloc((strlen(orig_name)+1)*sizeof(char));
				strcpy(pre_ctg[total].name,orig_name);
				++cul_id;
			}
			
			sscanf(line,">%s",orig_name);
			seq[0]='\0';
			len=0;
		}else{
			//strcat(seq,line);//effective?
			int single_len=strlen(line);
			line[single_len-1]='\0';
			strcpy(&seq[len],line);
			len+=single_len-1;
		}
	
	}
	if(len>overlaplen){
		boolean flag=upper_rev(seq,len);
		//fprintf(new_ctg,">%d length %d\n",cul_id,len);
		//fprintf(conver,"%s\t%d\t%d\n",orig_name,cul_id,len);
		//print_seq(new_ctg,seq,len);
		//fprintf(new_ctg,"%s\n",seq);
		char *one_seq=(char *)ckalloc((len+100)*sizeof(char));
		strcpy(one_seq,seq);
		if(flag==0){
			pre_ctg[++total].seq=one_seq;
			pre_ctg[total].bal=2;
			pre_ctg[total].len=len;
			pre_ctg[total].name=(char *)malloc(strlen(orig_name)*sizeof(char));
			strcpy(pre_ctg[total].name,orig_name);
			//pre_ctg[++total].bal=0;
			cul_id+=2;
		}else{
			pre_ctg[++total].seq=one_seq;
			pre_ctg[total].len=len;
			pre_ctg[total].bal=1;
			pre_ctg[total].name=(char *)malloc(strlen(orig_name)*sizeof(char));
			strcpy(pre_ctg[total].name,orig_name);
			++cul_id;
		}
		
	}
	fprintf(stderr,"All contigs loaded.\n");
	sprintf(file_name, "%s.updated.edge", graphfile);
	update=ckopen(file_name,"w");
	sprintf(file_name, "%s.ContigIndex", graphfile);
	index=ckopen(file_name,"w");
	fprintf(update,"EDGEs %d\n",cul_id);
	fprintf(index,"Edge_num %d %d\nindex\tlength\treverseComplement\n",cul_id,total);
	qsort(&pre_ctg[1],total,sizeof(IO_CTG),cmp_ctg);

	int i=1;
	cul_id=0;
	for(;i<=total;++i){
		if(pre_ctg[i].bal==2){
			len=pre_ctg[i].len;
			fprintf(new_ctg,">%d length %d\n",++cul_id,len);
			print_seq(new_ctg,pre_ctg[i].seq,len);
			fprintf(conver,"%s\t%d\t%d\n",pre_ctg[i].name,cul_id,len);
//			if(overlaplen<=31){
//				fprintf(update,">length %d,fffffffffff,fffffffffff,1,8\n",len);
//				fprintf(update,">length %d,fffffffffff,fffffffffff,-1,8\n",len);
//			}else{
				fprintf(update,">length %d,1,8\n",len);
				fprintf(update,">length %d,-1,8\n",len);
//			}
			fprintf(index,"%d\t%d\t1\n",cul_id++,len);
			
		}else{
			fprintf(new_ctg,">%d length %d\n",++cul_id,len);
			len=pre_ctg[i].len;
			print_seq(new_ctg,pre_ctg[i].seq,len);
			fprintf(conver,"%s\t%d\t%d\n",pre_ctg[i].name,cul_id,len);
			if(overlaplen<=31){
				fprintf(update,">length %d,fffffffffff,fffffffffff,0,8\n",len);
			}else{
				fprintf(update,">length %d,0,8\n",len);
			}
			fprintf(index,"%d\t%d\t0\n",cul_id,len);
		}
	}

	sprintf(file_name,"touch %s.Arc",graphfile);
	system(file_name);
	return 0;
}

//return value:0: in not equal its' rev_comp
//1: in equal its' rev_comp
boolean upper_rev(char *in,int in_len){
	int i,it_num;
	
	boolean ret_flag=1;
	it_num=in_len/2;
	
	for(i=0;i<it_num;++i){
		in[i]=toupper(in[i]);
		in[in_len-i-1]=toupper(in[in_len-i-1]);
		if(in[i]!=rev[in[in_len-i-1]]){
			ret_flag=0;
		}
	}
	return ret_flag;
}

void print_seq(FILE *out_file,char *sequence , int sequence_len){
	int it_num=sequence_len/100+1;
	int i;
	for(i=0;i<it_num;++i){
		char tmp;
		tmp=sequence[(i+1)*100];
		sequence[(i+1)*100]='\0';
		fprintf(out_file,"%s\n",&sequence[i*100]);
		sequence[(i+1)*100]=tmp;
	}
	

}

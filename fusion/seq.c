#include "stdinc.h"
#include "newhash.h"
#include "extfunc.h"
#include "extvab.h"

/*
put a insertSize in the grads array, 
if all grads have been entered and all the boundaris have been set, return 0
*/

void print_kmer(FILE *fp,Kmer kmer,char c)
{
	if(kmer)
		fprintf(fp,"%llx",kmer);
	else
		fprintf(fp,"0x0");
	fprintf(fp,"%c",c);

}

void printTightString(char *tightSeq,int len)
{
	int i;

	for(i=0;i<len;i++){
		printf("%c",int2base((int)getCharInTightString(tightSeq,i)));
		if((i+1)%100==0)
			printf("\n");
	}
	printf("\n");
}

static Kmer fastReverseComp(Kmer seq, char seq_size){
	seq ^= 0xAAAAAAAAAAAAAAAALLU;
	seq = ((seq & 0x3333333333333333LLU)<< 2) | ((seq & 0xCCCCCCCCCCCCCCCCLLU)>> 2);
	seq = ((seq & 0x0F0F0F0F0F0F0F0FLLU)<< 4) | ((seq & 0xF0F0F0F0F0F0F0F0LLU)>> 4);
	seq = ((seq & 0x00FF00FF00FF00FFLLU)<< 8) | ((seq & 0xFF00FF00FF00FF00LLU)>> 8);
	seq = ((seq & 0x0000FFFF0000FFFFLLU)<<16) | ((seq & 0xFFFF0000FFFF0000LLU)>>16);
	seq = ((seq & 0x00000000FFFFFFFFLLU)<<32) | ((seq & 0xFFFFFFFF00000000LLU)>>32);
	return seq >> (64 - (seq_size<<1));
}

Kmer reverseComplementVerbose(Kmer word,int overlap)
{
	return fastReverseComp(word,overlap);
	/*
	int index;
	Kmer revComp = 0;
	Kmer copy = word;
	unsigned char nucleotide;

	for (index = 0; index < overlap; index++) {
		nucleotide = copy & 3;
		revComp <<= 2;
		revComp += int_comp(nucleotide);//3 - nucleotide;
		copy >>= 2;
	}
	return revComp;
	*/
}

Kmer reverseComplement(Kmer word,int overlap)
{
	return fastReverseComp(word,overlap);
}

void writeChar2tightString(char nt,char *tightSeq,int pos)
{
	char *byte = tightSeq + pos/4;
	switch(pos%4){
	case 0:
		*byte &=63;
		*byte += nt << 6;
		return;	
	case 1:
		*byte &=207;
		*byte += nt << 4;
		return;	
	case 2:
		*byte &=243;
		*byte += nt << 2;
		return;	
	case 3:
		*byte &=252;
		*byte += nt;
		return;	
	
	}
}

char getCharInTightString(char *tightSeq,int pos)
{
	char *byte = tightSeq+pos/4;
	switch(pos%4){
	case 3:
		return (*byte & 3);
	case 2:
		return (*byte & 12) >> 2;
	case 1:
		return (*byte & 48) >> 4;
	case 0:
		return (*byte & 192) >> 6;
	}
	return 0;
}

// complement of sequence denoted 0, 1, 2, 3
void reverseComplementSeq(char *seq, int len,char *bal_seq)
{
	int i,index=0;
	
	if(len<1)
		return;

	for(i=len-1;i>=0;i--)
		bal_seq[index++] = int_comp(seq[i]);

	return;
}

// complement of sequence denoted 0, 1, 2, 3
char *compl_int_seq(char *seq, int len)
{
	char *bal_seq=NULL,c,bal_c;
	int i,index;

	if(len<1)
		return bal_seq;

	bal_seq = (char *)ckalloc(len*sizeof(char));
	index = 0;
	for(i=len-1;i>=0;i--){
		c = seq[i];
		if(c<4)
			bal_c = int_comp(c);//3-c;
		else 
			bal_c = c; 
		bal_seq[index++] = bal_c;

	}
	return bal_seq;
}

long long trans_seq(char *seq, int len)
{
        int     i;
        long long    res;

        res = 0;
        for(i = 0; i < len; i ++)       {
                res = res * 4 + seq[i];
        }

        return(res);
}

char *kmer2seq(Kmer word)
{
	int i;
	char *seq;
	Kmer charMask = 3;

	seq = (char *)ckalloc(overlaplen*sizeof(char));
	for(i=overlaplen-1;i>=0;i--){
		seq[i] = charMask&word;
		word >>= 2;
	}
	return seq;
}

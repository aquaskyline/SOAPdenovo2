#include "stdinc.h" 
#include "newhash.h"
#include "extfunc.h" 
#include "extvab.h" 

static unsigned char filter_array[8] = { (unsigned char) 1,((unsigned char) 1) << 1,((unsigned char) 1) << 2,((unsigned char) 1) << 3,((unsigned char) 1) << 4,((unsigned char) 1) << 5,((unsigned char) 1) << 6,((unsigned char) 1) << 7};


void link2next(NODE *node,char ch)
{	
	if(node->links & filter_array[(int)ch]) 
		node->linksB = node->linksB | filter_array[(int)ch]; 
	else
		node->links = node->links | filter_array[(int)ch]; 

}

unsigned char check_linkB2next(NODE *node,char ch)
{
	return filter_array[(int)ch]&node->linksB; 
}

unsigned char check_link2next(NODE *node,char ch)
{
	return filter_array[(int)ch]&node->links; 
}

void unlink2next(NODE *node,char ch)
{	
	node->links = node->links & (~filter_array[(int)ch]); 
}


void link2prev(NODE *node,char ch)
{	
	if(node->links & filter_array[ch+4]) 
		node->linksB = node->linksB | filter_array[ch+4]; 
	else
		node->links = node->links | filter_array[ch+4]; 
}

unsigned char check_linkB2prev(NODE *node,char ch)
{
	return filter_array[ch+4]&node->linksB; 
}

unsigned char check_link2prev(NODE *node,char ch)
{
	return filter_array[ch+4]&node->links; 
}

void unlink2prev(NODE *node,char ch)
{	
	node->links = node->links & (~filter_array[ch+4]); 
}

int count_link2next(NODE *node)
{
	int num = 0,i;
	unsigned char ch = node->links;

	for(i=0;i<4;i++){
		num += ch&0x01;
		ch >>= 1;
	}
	return num;
}

int count_link2nextB(NODE *node)
{
	int num = 0,i;
	unsigned char ch = node->linksB;

	for(i=0;i<4;i++){
		num += ch&0x01;
		ch >>= 1;
	}
	return num;
}

int count_link2prevB(NODE *node)
{
	int num = 0,i;
	unsigned char ch = node->linksB;

	ch >>= 4;
	for(i=0;i<4;i++){
		num += ch&0x01;
		ch >>= 1;
	}
	return num;
}

int count_link2prev(NODE *node)
{
	int num = 0,i;
	unsigned char ch = node->links;

	ch >>= 4;
	for(i=0;i<4;i++){
		num += ch&0x01;
		ch >>= 1;
	}
	return num;
}

Kmer KmerPlus(Kmer prev,char ch)
{
	Kmer word = prev;
	word <<= 2;
	word += ch;
	return word;
}
Kmer nextKmer(Kmer prev,char ch)
{
	Kmer word = prev;
	word <<= 2;
	word &= WORDFILTER;
	word += ch;
	return word;
}

Kmer prevKmer(Kmer next,char ch)
{
	Kmer word = next;
	word >>= 2;
	word += ((Kmer)ch) << 2*(overlaplen-1);
	return word;
}

char firstCharInKmer(Kmer kmer)
{
	return (char) (kmer >> 2*(overlaplen-1));// & 3;
}	


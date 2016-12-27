#include "stdinc.h"
#include "newhash.h"
#include "extfunc.h"
#include "extvab.h"

#define PUBLIC_FUNC
#define PROTECTED_FUNC

static const kmer_t empty_kmer = {0, 0, 0, 0, 0, 0, 1, 0,0};

static inline void update_kmer(kmer_t *mer, ubyte left, ubyte right){
	ubyte4 cov;

	if(left<4){
		cov = get_kmer_left_cov(*mer, left);
		if(cov < MAX_KMER_COV){
			set_kmer_left_cov(*mer, left, cov + 1);
		}
	}

	if(right<4){
		cov = get_kmer_right_cov(*mer, right);
		if(cov < MAX_KMER_COV){
			set_kmer_right_cov(*mer, right, cov + 1);
		}
	}
}

static inline void set_new_kmer(kmer_t *mer, ubyte8 seq, ubyte left, ubyte right){
	*mer = empty_kmer;
	set_kmer_seq(*mer, seq);
	if(left<4)
		set_kmer_left_cov(*mer, left, 1);
	if(right<4)
		set_kmer_right_cov(*mer, right, 1);
}


static inline int is_prime_kh(ubyte8 num){
	ubyte8 i, max;
	if(num < 4) return 1;
	if(num % 2 == 0) return 0;
	max = (ubyte8)sqrt((float)num);
	for(i=3;i<max;i+=2){ if(num % i == 0) return 0; }
	return 1;
}

static inline ubyte8 find_next_prime_kh(ubyte8 num){
	if(num % 2 == 0) num ++;
	while(1){ if(is_prime_kh(num)) return num; num += 2; }
}

PUBLIC_FUNC KmerSet* init_kmerset(ubyte8 init_size, float load_factor){
	KmerSet *set;
	if(init_size < 3) init_size = 3;
	else init_size = find_next_prime_kh(init_size);
	set = (KmerSet*)malloc(sizeof(KmerSet));
	set->size   = init_size;
	set->count  = 0;

	set->searchCnt = 0;
	set->foundCnt = 0;
	set->delCnt = 0;
	set->searchSpcSeedCnt = 0;
	set->getSpcSeedCnt = 0;
	set->levelGet[0] = 0;
	set->levelGet[1] = 0;
	set->levelGet[2] = 0;	
	
	set->max    = set->size * load_factor;
	if(load_factor <= 0) load_factor = 0.25f;
	else if(load_factor >= 1) load_factor = 0.75f;
	set->load_factor = load_factor;
	set->iter_ptr    = 0;
	set->array = calloc(set->size, sizeof(kmer_t));
	set->flags = malloc((set->size + 15)/16 * 4);
	memset(set->flags, 0x55, (set->size + 15) / 16 * 4);
	return set;
}

PROTECTED_FUNC static inline ubyte8 get_kmerset(KmerSet *set, ubyte8 seq){
	ubyte8 hc;
	hc = seq % set->size;
	while(1){
		if(is_kmer_entity_null(set->flags, hc)){
			return hc;
		} else {
			if(get_kmer_seq(set->array[hc]) == seq) return hc;
		}
		hc ++;
		if(hc == set->size) hc = 0;
	}
	return set->size;
}

PUBLIC_FUNC int search_kmerset(KmerSet *set, ubyte8 seq, kmer_t **rs){
	ubyte8 hc;
	hc = seq % set->size;
	while(1){
		if(is_kmer_entity_null(set->flags, hc)){
			return 0;
		} else {
			if(get_kmer_seq(set->array[hc]) == seq){
				*rs = set->array + hc;
				return 1;
			}
		}
		hc ++;
		if(hc == set->size) hc = 0;
	}
	return 0;
}

PUBLIC_FUNC static inline int exists_kmerset(KmerSet *set, ubyte8 seq){ 
	ubyte8 idx;
	idx = get_kmerset(set, seq);
	return !is_kmer_entity_null(set->flags, idx);
}

PROTECTED_FUNC static inline void encap_kmerset(KmerSet *set, ubyte8 num){
	ubyte4 *flags, *f;
	ubyte8 i, n, size, hc;
	kmer_t key, tmp;
	if(set->count + num <= set->max) return;
	n = set->size;
	do{
		if(n < 0xFFFFFFFU)
			n <<= 1;
		else
			n += 0xFFFFFFU;
		n = find_next_prime_kh(n); 
	} while(n * set->load_factor < set->count + num);

	set->array = realloc(set->array, n * sizeof(kmer_t));
	if(set->array == NULL){
		fprintf(stderr, "-- Out of memory --\n");
		abort();
	}
	flags = malloc((n+15)/16 * 4);
	memset(flags, 0x55, (n+15)/16 * 4);
	size = set->size;
	set->size = n;
	set->max = n * set->load_factor;
	f = set->flags;
	set->flags = flags;
	flags = f;
	for(i=0;i<size;i++){
		if(!exists_kmer_entity(flags, i)) continue;
		key = set->array[i];
		set_kmer_entity_del(flags, i);
		while(1){
			hc = get_kmer_seq(key) % set->size;
			while(!is_kmer_entity_null(set->flags, hc)){ hc ++; if(hc == set->size) hc = 0; }
			clear_kmer_entity_null(set->flags, hc);
			if(hc < size && exists_kmer_entity(flags, hc)){ 
				tmp = key;
				key = set->array[hc];
				set->array[hc] = tmp;
				set_kmer_entity_del(flags, hc);
			} else {
				set->array[hc] = key;
				break;
			}
		}
	}
	free(flags);
}

PUBLIC_FUNC int put_kmerset(KmerSet *set, ubyte8 seq, ubyte left, ubyte right, kmer_t **kmer_p){
	ubyte8 hc;
	encap_kmerset(set, 1);
	hc = seq % set->size;
	do{
		if(is_kmer_entity_null(set->flags, hc)){
			clear_kmer_entity_null(set->flags, hc);
			set_new_kmer(set->array + hc, seq, left, right);
			set->count ++;
			*kmer_p = set->array + hc;
			return 0;
		} else {
			if(get_kmer_seq(set->array[hc]) == seq){
				update_kmer(set->array + hc, left, right);
				set->array[hc].single = 0;
				*kmer_p = set->array + hc;
				return 1;
			}
		}
		hc ++;
		if(hc == set->size) hc = 0;
	} while(1);
	*kmer_p = NULL;
	return 0;
}

PUBLIC_FUNC byte8 count_kmerset(KmerSet *set){ return set->count; }

PUBLIC_FUNC static inline void reset_iter_kmerset(KmerSet *set){ set->iter_ptr = 0; }

PUBLIC_FUNC static inline ubyte8 iter_kmerset(KmerSet *set, kmer_t **rs){
	while(set->iter_ptr < set->size){
		if(!is_kmer_entity_null(set->flags, set->iter_ptr)){
			*rs = set->array + set->iter_ptr;
			set->iter_ptr ++;
			return 1;
		}
		set->iter_ptr ++;
	}
	return 0;
}

PUBLIC_FUNC void free_kmerset(KmerSet *set){
	free(set->array);
	free(set->flags);
	free(set);
}

PUBLIC_FUNC  void free_Sets(KmerSet **sets,int num){
	int i;
	for(i=0;i<num;i++)
		free_kmerset(sets[i]);
	free((void*)sets);
}

int count_branch2prev(kmer_t *node)
{
	int num = 0,i;

	for(i=0;i<4;i++){
		if(get_kmer_left_cov(*node,i)>0)
			num++;
	}
	return num;
}

int count_branch2next(kmer_t *node)
{
	int num = 0,i;

	for(i=0;i<4;i++){
		if(get_kmer_right_cov(*node,i)>0)
			num++;
	}
	return num;
}

void dislink2prevUncertain(kmer_t *node,char ch,boolean smaller)
{
	if(smaller)
		set_kmer_left_cov(*node,ch,0);
	else
		set_kmer_right_cov(*node,int_comp(ch),0);
	
}

void dislink2nextUncertain(kmer_t *node,char ch,boolean smaller)
{
	if(smaller)
		set_kmer_right_cov(*node,ch,0);
	else
		set_kmer_left_cov(*node,int_comp(ch),0);
}
	





////////////////// functions for spaced seed Kmer hash

static const spcKmer empty_spckmer = {0, NULL, 1};

static inline int update_spckmer(spcKmer *mer, ubyte2 s_bases, kmer_t *node){
//	if(mer->start == NULL)
//		fprintf(stderr, "start err at:\t%llu\n",mer->seq);
		
	spcBase *tmpBase=mer->start;

	spcBase *newSpcBase;
	newSpcBase = (spcBase*)malloc(sizeof(spcBase));
	newSpcBase->spaced_bases = s_bases;
//	newSpcBase->edgeID = edgeID;
	newSpcBase->large_kmer = node;
	newSpcBase->next = tmpBase->next;
	tmpBase->next = newSpcBase;
	
	mer->spaced_base_num++;
	
//	mvnv(0,"update %llu :\t%hu\tnum: %u\n", mer->seq, tmpBase->next->spaced_bases, mer->spaced_base_num);
	return 0;
}

static inline void set_new_spckmer(spcKmer *mer, Kmer spc_kmer, ubyte2 s_bases, kmer_t *node){
	*mer = empty_spckmer;
	set_kmer_seq(*mer, spc_kmer);
	
	spcBase *newSpcBase;
	newSpcBase = (spcBase*)malloc(sizeof(spcBase));
	newSpcBase->spaced_bases = s_bases;
//	newSpcBase->repeat = 0;
//	newSpcBase->edgeID = edgeID;
	newSpcBase->large_kmer = node;
	newSpcBase->next = NULL;
	
	mer->start = newSpcBase;
	
//	mvnv(0,"new %llu :\t%hu\n", mer->seq, mer->start->spaced_bases)
	
}

PUBLIC_FUNC spcKmerSet* init_spckmerset(ubyte8 init_size, float load_factor){
	spcKmerSet *set;
	if(init_size < 3) init_size = 3;
	else init_size = find_next_prime_kh(init_size);

	set = (spcKmerSet*)malloc(sizeof(spcKmerSet));
	set->size   = init_size;
	set->count  = 0;
	set->max    = set->size * load_factor;
	if(load_factor <= 0) load_factor = 0.25f;
	else if(load_factor >= 1) load_factor = 0.75f;
	set->load_factor = load_factor;
	//set->iter_ptr    = 0;
	set->array = calloc(set->size, sizeof(spcKmer));
	set->flags = malloc((set->size + 15)/16 * 4);
	memset(set->flags, 0x55, (set->size + 15) / 16 * 4);
	return set;
}

PUBLIC_FUNC int search_spckmerset(spcKmerSet *set, ubyte8 seq, spcKmer **rs){
	ubyte8 hc;
	hc = seq % set->size;
	while(1){
		if(is_kmer_entity_null(set->flags, hc)){
			return 0;
		} else {
			if(get_kmer_seq(set->array[hc]) == seq){
				*rs = set->array + hc;
				return 1;
			}
		}
		hc ++;
		if(hc == set->size) hc = 0;
	}
	return 0;
}

PROTECTED_FUNC static inline void encap_spckmerset(spcKmerSet *set, ubyte8 num){
	ubyte4 *flags, *f;
	ubyte8 i, n, size, hc;
	spcKmer key, tmp;
	if(set->count + num <= set->max) return;

	n = set->size;
	do{
		if(n < 0xFFFFFFFU)
			n <<= 1;
		else
			n += 0xFFFFFFU;
		n = find_next_prime_kh(n); 
	} while(n * set->load_factor < set->count + num);

	set->array = realloc(set->array, n * sizeof(spcKmer));
	if(set->array == NULL){
		fprintf(stderr, "-- Out of memory --\n");
		abort();
	}
	flags = malloc((n+15)/16 * 4);
	memset(flags, 0x55, (n+15)/16 * 4);
	size = set->size;
	set->size = n;
	set->max = n * set->load_factor;
	f = set->flags;
	set->flags = flags;
	flags = f;
	for(i=0;i<size;i++){
		if(!exists_kmer_entity(flags, i)) continue;
		key = set->array[i];
		set_kmer_entity_del(flags, i);
		while(1){
			hc = get_kmer_seq(key) % set->size;
			while(!is_kmer_entity_null(set->flags, hc)){ hc ++; if(hc == set->size) hc = 0; }
			clear_kmer_entity_null(set->flags, hc);
			if(hc < size && exists_kmer_entity(flags, hc)){ 
				tmp = key;
				key = set->array[hc];
				set->array[hc] = tmp;
				set_kmer_entity_del(flags, hc);
			} else {
				set->array[hc] = key;
				break;
			}
		}
	}
	free(flags);
}

PUBLIC_FUNC int put_spckmerset(spcKmerSet *set, Kmer spc_kmer, ubyte2 spaced_bases, kmer_t *node){
	ubyte8 hc;
	encap_spckmerset(set, 1);
	hc = spc_kmer % set->size;
	do{
		if(is_kmer_entity_null(set->flags, hc)){		//new! repeat_flag==0
			clear_kmer_entity_null(set->flags, hc);
			set_new_spckmer(set->array + hc, spc_kmer, spaced_bases, node);
			set->count ++;
			return 0;
		} else {
			if(get_kmer_seq(set->array[hc]) == spc_kmer){	//exists! repeat_flag==1 or 0
				return update_spckmer(set->array + hc, spaced_bases, node);
			}
		}
		hc ++;
		if(hc == set->size) hc = 0;
	} while(1);
	return 3;
}

PUBLIC_FUNC void buildSpcKmerSet(KmerSet *set, spcKmerSet *spaced_kset)
{
	boolean spcFlag;
	Kmer buff_kmer, spc_kmer;
	ubyte2 spc_bases;
	
	ubyte8 i=0,j=0;
	for(i=0;i<set->size;i++)
	{
		if(is_kmer_entity_null(set->flags, i))
			continue;
		else
		{
//			kmer_t **kmer_p;
//			*kmer_p = set->array+i;
			if(set->array[i].deleted != 1)	//kmer not repeat
			{
				//spaced seed: 18 of 25, build masker and use >>,&,| for each part, only assign once
				//	1 1111 1010 1100 1111 1101 0110			!!!OLD!!!
				//	1 1111 1111 1111 1010 1100 1000			!!!NEW!!!
				//	11 11111111 11111111 11111111 11001100 11110000 11000000			!!!NEW!!!				
				
				buff_kmer = get_kmer_seq(set->array[i]);
				
				spc_kmer = ((buff_kmer>>14)&0xFFFFFFF00) | ((buff_kmer>>12)&0xC0) | ((buff_kmer>>10)&0x3C) | ((buff_kmer>>6)&0x3);
				//0xFFFFFFF00 = 1111 11111111 11111111 11111111 00000000
				//			 0xC0 = 0000 00000000 00000000 00000000 11000000
				//			 0x3C = 0000 00000000 00000000 00000000 00111100
				//				0x3 = 0000 00000000 00000000 00000000 00000011
				
				spc_bases = ((buff_kmer>>8)&0x3000) | ((buff_kmer>>6)&0xC00) | ((buff_kmer>>2)&0x3C0) | (buff_kmer&0x3F);
				//		 0x3000 = 110000 00000000
				//			0xC00 = 001100 00000000
				//			0x3C0 = 000011 11000000
				//			 0x3F = 000000 00111111

				//build the 18mer and the spaced bases(7mer), put them in the spaced_kmer hash				
				spcFlag = put_spckmerset(spaced_kset, spc_kmer, spc_bases, set->array+i);
				if(spcFlag!=0)
					fprintf(stderr, "flag error: %c\tkmer exists: %llu %hu\n", spcFlag, spc_kmer, spc_bases);
//				if((++j)%100000==0)
//					fprintf(stderr,"--- %lluth spaced Kmer built\n",j);		
			}
			
		}
	}
	//fprintf(stderr,"--- total %llu spaced Kmer built in a KmerSet\n",j);	
}

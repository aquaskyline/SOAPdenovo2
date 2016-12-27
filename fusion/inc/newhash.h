#ifndef __NEW_HASH_RJ
#define __NEW_HASH_RJ

#ifndef K_LOAD_FACTOR
#define K_LOAD_FACTOR 0.75
#endif

#define MAX_KMER_COV 63
#define EDGE_BIT_SIZE 6
#define EDGE_XOR_MASK 0x3FU
#define LINKS_BITS 0x00FFFFFFU

#define get_kmer_seq(mer) ((mer).seq)
#define set_kmer_seq(mer, val) ((mer).seq = val)

#define get_kmer_left_cov(mer, idx) (((mer).l_links>>((idx)*EDGE_BIT_SIZE))&EDGE_XOR_MASK)
#define set_kmer_left_cov(mer, idx, val) ((mer).l_links = ((mer).l_links&(~(EDGE_XOR_MASK<<((idx)*EDGE_BIT_SIZE)))) | (((val)&EDGE_XOR_MASK)<<((idx)*EDGE_BIT_SIZE)) )
#define get_kmer_left_covs(mer) (get_kmer_left_cov(mer, 0) + get_kmer_left_cov(mer, 1) + get_kmer_left_cov(mer, 2) + get_kmer_left_cov(mer, 3))

#define get_kmer_right_cov(mer, idx) (((mer).r_links>>((idx)*EDGE_BIT_SIZE))&EDGE_XOR_MASK)
#define set_kmer_right_cov(mer, idx, val) ((mer).r_links = ((mer).r_links&(~(EDGE_XOR_MASK<<((idx)*EDGE_BIT_SIZE)))) | (((val)&EDGE_XOR_MASK)<<((idx)*EDGE_BIT_SIZE)) )
#define get_kmer_right_covs(mer) (get_kmer_right_cov(mer, 0) + get_kmer_right_cov(mer, 1) + get_kmer_right_cov(mer, 2) + get_kmer_right_cov(mer, 3))


#define is_kmer_entity_null(flags, idx)    ((flags)[(idx)>>4]>>(((idx)&0x0f)<<1)&0x01)
#define is_kmer_entity_del(flags, idx)     ((flags)[(idx)>>4]>>(((idx)&0x0f)<<1)&0x02)
#define set_kmer_entity_null(flags, idx)   ((flags)[(idx)>>4] |= (0x01u<<(((idx)&0x0f)<<1)))
#define set_kmer_entity_del(flags, idx)    ((flags)[(idx)>>4] |= (0x02u<<(((idx)&0x0f)<<1)))
#define clear_kmer_entity_null(flags, idx) ((flags)[(idx)>>4] &= ~(0x01u<<(((idx)&0x0f)<<1)))
#define clear_kmer_entity_del(flags, idx)  ((flags)[(idx)>>4] &= ~(0x02u<<(((idx)&0x0f)<<1)))
#define exists_kmer_entity(flags, idx)     (!((flags)[(idx)>>4]>>(((idx)&0x0f)<<1)&0x03))


typedef struct kmer_st
{
	Kmer seq;
	ubyte4 l_links;                     // sever as edgeID since make_edge
	ubyte4 r_links:4*EDGE_BIT_SIZE;
	ubyte4 linear:1;
	ubyte4 deleted:1;
	ubyte4 checked:1;
	ubyte4 single:1;
	ubyte4 twin:2;
	ubyte4 inEdge:2;
} kmer_t;

typedef struct kmerSet_st
{
	kmer_t *array;
	ubyte4 *flags;
	ubyte8 size;
	ubyte8 count;
	ubyte8 max;
	double load_factor;
	ubyte8 iter_ptr;
	
	ubyte8 searchCnt;
	ubyte8 foundCnt;
	ubyte8 delCnt;
	ubyte8 searchSpcSeedCnt;
	ubyte8 getSpcSeedCnt;
	ubyte8 levelGet[3];
	
} KmerSet;

typedef struct kmer_pt
{
	kmer_t *node;
	Kmer kmer;
	boolean isSmaller;
	struct kmer_pt *next;
}KMER_PT;

//////////////////////////////////////////////////////////////// spaced seed

typedef struct spaced_base
{
	ubyte2 spaced_bases:14;
	//ubyte2 repeat:1;
	//ubyte4 edgeID;
	kmer_t *large_kmer;
	struct spaced_base *next;
}spcBase;

typedef struct spaced_kmer
{
	Kmer seq;
	struct spaced_base *start;
	ubyte4 spaced_base_num;
}spcKmer;

typedef struct spaced_kmer_set
{
	spcKmer *array;
	ubyte4 *flags;
	ubyte8 size;
	ubyte8 count;
	ubyte8 max;
	double load_factor;
} spcKmerSet;

extern spcKmerSet* init_spckmerset(ubyte8 init_size, float load_factor);
extern void buildSpcKmerSet(KmerSet *set, spcKmerSet *spaced_kset);
extern int search_spckmerset(spcKmerSet *set, ubyte8 seq, spcKmer **rs);
extern int put_spckmerset(spcKmerSet *set, Kmer spc_kmer, ubyte2 spaced_bases, kmer_t *node);

//////////////////////////////////////////////////////////////// spaced seed END

extern KmerSet* init_kmerset(ubyte8 init_size, float load_factor);
extern int search_kmerset(KmerSet *set, ubyte8 seq, kmer_t **rs);
extern int put_kmerset(KmerSet *set, ubyte8 seq, ubyte left, ubyte right,kmer_t **kmer_p);
extern byte8 count_kmerset(KmerSet *set);
extern void free_Sets(KmerSet **KmerSets,int num);
extern void free_kmerset(KmerSet *set);
extern void dislink2nextUncertain(kmer_t *node,char ch,boolean smaller);
extern void dislink2prevUncertain(kmer_t *node,char ch,boolean smaller);

extern int count_branch2prev(kmer_t *node);
extern int count_branch2next(kmer_t *node);
extern char firstCharInKmer(Kmer kmer);

#endif

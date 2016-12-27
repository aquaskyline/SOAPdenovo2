#ifndef __STACK__
#define __STACK__
 
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

typedef struct block_starter
{
	struct block_starter *prev;
	struct block_starter *next;
}BLOCK_STARTER;

typedef struct stack
{
	BLOCK_STARTER *block_list;
	int index_in_block;
	int items_per_block;
	int item_c;
	size_t item_size;
	BLOCK_STARTER *block_backup;
	int index_backup;
	int item_c_backup;
}STACK;

void stackBackup(STACK *astack);
void stackRecover(STACK *astack);
void *stackPush(STACK *astack);
void *stackPop(STACK *astack);
void freeStack(STACK *astack);
void emptyStack(STACK *astack);
STACK *createStack(int num_items,size_t unit_size);


#endif

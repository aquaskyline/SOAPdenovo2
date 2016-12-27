#include "stack.h"

STACK *createStack(int num_items,size_t unit_size)
{
	STACK *newStack = (STACK *)malloc(1*sizeof(STACK));

	newStack->block_list = NULL;
	newStack->items_per_block = num_items;
	newStack->item_size = unit_size;
	newStack->item_c = 0;
	return newStack; 
}

void emptyStack(STACK *astack)
{
	BLOCK_STARTER *block;
	if(!astack||!astack->block_list)
		return;

	block = astack->block_list;
	if(block->next)
		block = block->next;

	astack->block_list = block;
	astack->item_c = 0;
	astack->index_in_block = 0;
}

void freeStack(STACK *astack)
{
	BLOCK_STARTER *ite_block,*temp_block;

	if(!astack)
		return;

	ite_block = astack->block_list;
	if(ite_block){
		while(ite_block->next)
			ite_block = ite_block->next;
	}
	while(ite_block){
		temp_block = ite_block;
		ite_block = ite_block->prev;
		free((void *)temp_block);
	}

	free((void *)astack);
}

void stackBackup(STACK *astack)
{
	astack->block_backup = astack->block_list;
	astack->index_backup = astack->index_in_block;
	astack->item_c_backup = astack->item_c;
}

void stackRecover(STACK *astack)
{
	astack->block_list = astack->block_backup;
	astack->index_in_block = astack->index_backup;
	astack->item_c = astack->item_c_backup;
}

void *stackPop(STACK *astack)
{
	BLOCK_STARTER *block;

	if(!astack||!astack->block_list||!astack->item_c)
		return NULL;
	
	astack->item_c--;
	block = astack->block_list;
	if(astack->index_in_block==1){
		if(block->next){
			astack->block_list = block->next;
			astack->index_in_block = astack->items_per_block;
		}else{
			astack->index_in_block = 0;
			astack->item_c = 0;
		}
		return (void *)((void *)block+sizeof(BLOCK_STARTER));

	}
	return (void *)((void *)block+sizeof(BLOCK_STARTER)+astack->item_size*(--astack->index_in_block));
}

void *stackPush(STACK *astack)
{
	BLOCK_STARTER *block;

	if(!astack)
		return NULL;

	astack->item_c++;
	if(!astack->block_list||(astack->index_in_block==astack->items_per_block&&!astack->block_list->prev)){
		block = malloc(sizeof(BLOCK_STARTER)+astack->items_per_block*astack->item_size);
		block->prev = NULL;
		if(astack->block_list)
			astack->block_list->prev = block;
		block->next = astack->block_list;
		astack->block_list = block;
		astack->index_in_block = 1;
		return (void *)((void *)block+sizeof(BLOCK_STARTER));
	}else if(astack->index_in_block==astack->items_per_block&&astack->block_list->prev){
		astack->block_list = astack->block_list->prev;
		astack->index_in_block = 1;
		return (void *)((void *)astack->block_list+sizeof(BLOCK_STARTER));
	}

	block = astack->block_list;
	return (void *)((void *)block+sizeof(BLOCK_STARTER)+astack->item_size*astack->index_in_block++);
	
}

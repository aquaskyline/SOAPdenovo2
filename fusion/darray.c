#include "darray.h"
#include "check.h"

DARRAY *createDarray(int num_items,size_t unit_size)
{
	DARRAY *newDarray = (DARRAY *)malloc(1*sizeof(DARRAY));

	newDarray->array_size = num_items;
	newDarray->item_size = unit_size;
	newDarray->item_c = 0;
	newDarray->array = (void *)ckalloc(num_items*unit_size);
	return newDarray; 
}

void *darrayPut(DARRAY *darray,long long index)
{
	int i=2;
	if(index+1>darray->item_c)
		darray->item_c = index + 1;
	if(index<darray->array_size)
		return darray->array + darray->item_size*index;
	while(index>i*darray->array_size)
		i++;
	
	darray->array = (void *)ckrealloc(darray->array,i*darray->array_size*darray->item_size
				,darray->array_size*darray->item_size);
	darray->array_size *=i;
	return (void *)((void *)darray->array + darray->item_size*index);
}

void *darrayGet(DARRAY *darray, long long index)
{
	if(index<darray->array_size)
		return (void *)((void *)darray->array + darray->item_size*index);
	printf("array read index %lld out of range %lld\n",index,darray->array_size);
	return NULL;
}


void emptyDarray(DARRAY *darray)
{
	darray->item_c = 0;
}

void freeDarray(DARRAY *darray)
{

	if(!darray)
		return;

	if(darray->array)
		free((void *)darray->array);

	free((void *)darray);
}


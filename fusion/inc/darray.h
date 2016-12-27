#ifndef __DARRAY__
#define __DARRAY__

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

typedef struct dynamic_array 
{
	void *array;
	long long array_size;
	size_t item_size;
	long long item_c;
}DARRAY;

void *darrayPut(DARRAY *darray,long long index);
void *darrayGet(DARRAY *darray,long long index);
DARRAY *createDarray(int num_items,size_t unit_size);
void freeDarray(DARRAY *darray);
void emptyDarray(DARRAY *darray);

#endif


/***************************************************************************
 * Title:          check.c
 * Author:         Haixu Tang
 * Created:        Jun. 2002
 * Last modified:  May. 2004
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

/* ckopen - open file; check for success */

#include <stdinc.h>
//#include <extfunc.h>

void *ckalloc(unsigned long long amount);
FILE *ckopen(char *name, char *mode);

FILE *ckopen(char *name, char *mode)
{
	FILE *fp;

	if ((fp = fopen(name, mode)) == NULL)	{
		printf("Cannot open file %s.\n", name);
		exit(-1);
	}
	return(fp);
}


/* ckalloc - allocate space; check for success */

void *ckalloc(unsigned long long amount)
{
	void *p;

	if ((p = (void *) calloc( 1, (unsigned long long) amount)) == NULL && amount != 0)	{
		printf("not enought memory");
		fflush(stdout);
                exit(-1);
	}
	return(p);
}


/* reallocate memory */
void *ckrealloc(void *p, size_t new_size, size_t old_size)
{
  void *q;

  q = realloc((void *) p, new_size);
  if (new_size == 0 || q != (void *) 0)
    return q;

  /* manually reallocate space */
  q = ckalloc(new_size);

  /* move old memory to new space */
  bcopy(p, q, old_size);
  free(p);

  return q;
}

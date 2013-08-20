/*
 * readInterval.c
 *
 * Copyright (c) 2008-2012 BGI-Shenzhen <soap at genomics dot org dot cn>.
 *
 * This file is part of SOAPdenovo.
 *
 * SOAPdenovo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SOAPdenovo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SOAPdenovo.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "stdinc.h"
#include "newhash.h"
#include "kmerhash.h"
#include "extfunc.h"
#include "extvab.h"

#define RVBLOCKSIZE 1000

void destroyReadIntervMem ()
{
	freeMem_manager ( rv_mem_manager );
	rv_mem_manager = NULL;
}

READINTERVAL * allocateRV ( int readid, int edgeid )
{
	READINTERVAL * newRV;
	newRV = ( READINTERVAL * ) getItem ( rv_mem_manager );
	newRV->readid = readid;
	newRV->edgeid = edgeid;
	newRV->nextInRead = NULL;
	newRV->prevInRead = NULL;
	newRV->nextOnEdge = NULL;
	newRV->prevOnEdge = NULL;
	return newRV;
}

void dismissRV ( READINTERVAL * rv )
{
	returnItem ( rv_mem_manager, rv );
}

/*************************************************
Function:
    createRVmemo
Description:
    Alloc the memory for rv_mem_manager.
Input:
    None.
Output:
    None.
Return:
    None.
*************************************************/
void createRVmemo ()
{
	if ( !rv_mem_manager )
	{
		rv_mem_manager = createMem_manager ( RVBLOCKSIZE, sizeof ( READINTERVAL ) );
	}
	else
	{
		fprintf ( stderr, "Warning from createRVmemo: rv_mem_manager is an active pointer.\n" );
	}
}

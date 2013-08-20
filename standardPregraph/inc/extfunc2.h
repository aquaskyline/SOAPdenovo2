/*
 * inc/extfunc2.h
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

#ifndef _MEM_MANAGER
#define _MEM_MANAGER
extern MEM_MANAGER * createMem_manager ( int num_items, size_t unit_size );
extern void * getItem ( MEM_MANAGER * mem_Manager );
extern void returnItem ( MEM_MANAGER * mem_Manager, void * );
extern void freeMem_manager ( MEM_MANAGER * mem_Manager );
#endif

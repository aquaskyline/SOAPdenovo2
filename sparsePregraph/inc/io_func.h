/*
 * inc/io_func.h
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


#ifndef _IO_FUNC_H
#define _IO_FUNC_H

#include "stdinc.h"
#include "global.h"

struct io_para_main
{
	int read_buf_sz;
	vector<string> *in_filenames_vt;

};

void sendIOWorkSignal();
void * run_io_thread_main ( void * arg );
void filter_N ( string & seq_s, int & bad_flag );
void read_lib ( vector<string> &filenames, char * lib_file );
void * open_file_robust ( const char * filetype, const char * path, const char * mode );

#endif











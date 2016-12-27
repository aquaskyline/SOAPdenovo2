/*
 *      Filename: general.h
 *
 *
 *      Description:
 *         Basic functions
 *
 *  	Created on: Feb 8, 2010
 *      Author: Ruibang Luo, BGI
 *
 *     	History:
 *         1.
 */

#pragma once
#ifndef GENERAL_H_AQUA_
#define GENERAL_H_AQUA_

#include<unistd.h>

//Useful Variables*************************************************************
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#define FN_SIZE 2048
//*****************************************************************************

//Types************************************************************************
typedef unsigned int       uint;
typedef unsigned char      uchar;
typedef unsigned short     ushort;
typedef unsigned long      ulong;
typedef unsigned long long ulonglong;

typedef unsigned char      BYTE;
typedef unsigned short     WORD;
typedef unsigned int       DWORD;

typedef unsigned char      u8_t;
typedef unsigned short     u16_t;
typedef unsigned int       u32_t;
typedef unsigned long long u64_t;

typedef char * chptr;

//*****************************************************************************

//Debugging********************************************************************
//Verbose system
//Verbosity should seperated into 4 levels: 0, 1, 2, 3
#define VERBOSITY_BOTTOM 0
#define VERBOSITY_TOP 4
int ModifyVerbosity(const int);
#define verboseBufSize 16384

#define ModVerboseStrAndVerbose(level, ...) \
			{\
				if(verbosity >> level)\
				{\
					snprintf(verboseStr, verboseBufSize, ##__VA_ARGS__);\
					fprintf(stderr,"[%s]:%s\n",__FUNCTION__,verboseStr);\
				}\
			}
#define mvnv(level, ...) ModVerboseStrAndVerbose(level, ##__VA_ARGS__)
#define die(...) \
		{\
			ModVerboseStrAndVerbose(0, ##__VA_ARGS__);\
			fprintf(stderr,"Program terminated.\n");\
			exit(EXIT_FAILURE);\
		}
#define sigdie(sig, ...) \
		{\
			ModVerboseStrAndVerbose(0, ##__VA_ARGS__);\
			fprintf(stderr,"Program terminated.\n");\
			exit(sig);\
		}
#define perrdie(...) \
		{\
			ModVerboseStrAndVerbose(0, ##__VA_ARGS__);\
			perror("");\
			fprintf(stderr,"Program terminated.\n");\
			exit(EXIT_FAILURE);\
		}
#define mk \
{\
	fprintf(stderr, "DBG Marker @ %s:%d\n", __FUNCTION__, __LINE__);\
}

#endif

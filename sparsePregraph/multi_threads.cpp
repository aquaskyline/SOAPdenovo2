/*
 * multi_threads.cpp
 *
 * Copyright (c) 2008-2016 Ruibang Luo <aquaskyline.com>.
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

#include "multi_threads.h"
#include "global.h"
#include "build_graph.h"
#include "stdinc.h"

#include <unistd.h>

#include "build_preArc.h"

void creatThrds ( pthread_t *threads, PARAMETER *paras )
{
  unsigned char i;
  int temp;

  for ( i = 0; i < thrd_num_s; i++ )
    {
      if ( ( temp = pthread_create ( &threads[i], NULL, threadRoutine, & ( paras[i] ) ) ) != 0 )
        {
          fprintf ( stderr, "ERROR: create threads failed.\n" );
          exit ( 1 );
        }
    }

  //fprintf(stderr,"%d work threads created.\n",thrd_num_s);
  fprintf ( stderr, "%d work threads initialized.\n", thrd_num_s );
}



void *threadRoutine ( void *para )
{
  PARAMETER *prm;
  int i;
  unsigned char id;
  prm = ( PARAMETER * ) para;
  id = prm->threadID;
  struct hashtable2 *ht = prm->ht;
  int K_size = prm->K_size;
  int gap = prm->gap;

  while ( 1 )
    {
      if ( * ( prm->selfSignal ) == 3 )
        {
          * ( prm->selfSignal ) = 0;
          break;
        }
      else if ( * ( prm->selfSignal ) == 10 )
        {
          run_process_threaded ( ht, locks, K_size, gap, read_num, thrd_num_s, prm->threadID, 1 );
          * ( prm->selfSignal ) = 0;
        }
      else if ( * ( prm->selfSignal ) == 11 )
        {
          run_process_threaded ( ht, locks, K_size, gap, read_num, thrd_num_s, prm->threadID, 2 );
          * ( prm->selfSignal ) = 0;
        }
      else if ( * ( prm->selfSignal ) == 12 )
        {
          for ( int i = prm->threadID ; i < read_num; i += thrd_num_s )
            {
              int bad_flag = 0;
              filter_N ( seq_t[i], bad_flag );

              if ( bad_flag )
                {
                  seq_t[i].clear();
                  continue;
                }

              process_1read_preArc ( prm->preArcs, locks, prm->threadID, prm->v_ht, K_size, prm->cut_off_len, seq_t[i].c_str() );
            }

          * ( prm->selfSignal ) = 0;
        }

      usleep ( 1 );
    }
}

void thread_wait ( pthread_t *threads )
{
  int i;

  for ( i = 0; i < thrd_num_s; i++ )
    if ( threads[i] != 0 )
      {
        pthread_join ( threads[i], NULL );
      }
}

void sendWorkSignal ( unsigned char SIG, unsigned char *thrdSignals )
{
  int t;

  for ( t = 0; t < thrd_num_s; t++ )
    {
      thrdSignals[t + 1] = SIG;
    }

  while ( 1 )
    {
      usleep ( 10 );

      for ( t = 0; t < thrd_num_s; t++ )
        if ( thrdSignals[t + 1] )
          {
            break;
          }

      if ( t == thrd_num_s )
        {
          break;
        }
    }
}






















/*
 * pregraph_sparse.cpp
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

#include "global.h"
#include "stdinc.h"
#include "core.h"

#include <unistd.h>

#include "multi_threads.h"
#include "build_graph.h"
#include "build_edge.h"
#include "build_preArc.h"
#include "io_func.h"
#include "seq_util.h"
#include "convert_soapdenovo.h"


static int LOAD_GRAPH = 0, BUILD_DBG = 1, BUILD_EDGES = 1, BUILD_PREARCS = 1;
//static    int run_mode=0;

extern "C" int call_pregraph_sparse ( int argc, char **argv );
static void initenv ( int argc, char **argv );
static void parse_args ( vector<string> &in_filenames_vt );
static void display_pregraph_usage();

/*
int main ( int argc, char ** argv )
{
	fprintf ( stderr, "\nVersion 1.0.3: released on July 13th, 2012\nCompile %s\t%s\n\n", __DATE__, __TIME__ );
	call_pregraph_sparse ( argc, argv );
}*/



/*************************************************
Function:
    call_pregraph_sparse
Description:
    The entry function for calling pregraph_sparse, its processes include
    1. Parses the args form command.
    2. Builds the sparse-kmer graph.
    3. Removes low coverage sparse-kmers and kmer-edges (the connections between edges).
    4. Removes tips.
    5. Builds edges by compacting linear kmer nodes.
    6. Builds preArc (the connections between edges) by maping reads to edges.
Input:
    @see  display_pregraph_usage()
Output:
    files:
    *.ht_idx
    *.ht_content
    *.kmerFreq
    *.sparse.edge
    *.edge.gz
    *.preArc
Return:
    None.
*************************************************/
extern "C" int call_pregraph_sparse ( int argc, char **argv )
{
  time_t all_beg_time, all_end_time;
  time ( &all_beg_time );
  initenv ( argc, argv );
  vector<string> in_filenames_vt;
  parse_args ( in_filenames_vt );
  struct hashtable2 ht2;
  size_t hashTableSZ = 1000000;
  time_t beg_time, read_time;
  size_t bucket_count = 0, edge_cnt = 0;

  if ( ( !LOAD_GRAPH ) && BUILD_DBG )
    {
      int round = 1;

      for ( round = 1; round <= 2; round++ )
        {
          //fprintf(stderr,"Building the sparse de Brujin graph, round: %d\n",round);
          fprintf ( stderr, "Start to build the sparse de Brujin graph, round: %d\n", round );

          if ( round == 1 )
            {
              time ( &beg_time );
              edge_cnt = 0;
              uint64_t TotalSamplings = 0;

              //initialize the hashtable size
              if ( GenomeSize == 0 )
                {
                  fprintf ( stderr, "Error! Genome size not given.\n" );
                  return -1;
                }

              hashTableSZ = ( size_t ) GenomeSize / max ( gap - 1, 5 );
              int read_buf_sz = 102400 * thrd_num_s;
              Init_HT2 ( &ht2, hashTableSZ );
              //create main io thread
              read_buf0 = new string[read_buf_sz];
              read_buf1 = new string[read_buf_sz];
              io_stat0 = 1; //must be one, if io_stat0 equals 0, the io threads will work immediately
              io_stat1 = 1;
              io_ready = 0;
              io_para_main io_para_mains;
              io_para_mains.read_buf_sz = read_buf_sz;
              io_para_mains.in_filenames_vt = &in_filenames_vt;
              pthread_t io_thread;
              int temp;

              // fprintf(stderr,"Creating main io thread ...\n");
              if ( ( temp = pthread_create ( &io_thread, NULL, run_io_thread_main, &io_para_mains ) ) != 0 )
                {
                  fprintf ( stderr, "ERROR: failed creating main io thread.\n" );
                  exit ( -1 );
                }

              fprintf ( stderr, "1 io thread initialized.\n" );
              //create work threads  for round 1
              pthread_t threads[thrd_num_s];
              unsigned char thrdSignal[thrd_num_s + 1];
              PARAMETER paras[thrd_num_s];
              locks = ( pthread_spinlock_t * ) calloc ( ht2.ht_sz, sizeof ( pthread_spinlock_t ) );

              //initialize the locks  unlock
              for ( size_t i = 0; i < ht2.ht_sz; ++i )
                {
                  locks[i] = 1;
                }

              //create threads
              //fprintf(stderr,"Creating work threads ...\n");
              bucket_count_total = ( size_t * ) calloc ( thrd_num_s, sizeof ( size_t ) );
              edge_cnt_total = ( size_t * ) calloc ( thrd_num_s, sizeof ( size_t ) );

              for ( int k = 0; k < thrd_num_s; k++ )
                {
                  thrdSignal[k + 1] = 0;
                  paras[k].threadID = k;
                  paras[k].mainSignal = &thrdSignal[0];
                  paras[k].selfSignal = &thrdSignal[k + 1];
                  paras[k].ht = &ht2;
                  paras[k].K_size = K_size;
                  paras[k].gap = gap;
                }

              creatThrds ( threads, paras );
              thrdSignal[0] = 0;
              //begin to work
              size_t processed_reads = 0;

              while ( 1 )
                {
                  sendIOWorkSignal();

                  while ( io_ready == 0 )
                    {
                      usleep ( 1 );
                    }

                  if ( io_ready )
                    {
                      processed_reads += read_num;
                      sendWorkSignal ( 10, thrdSignal );

                      for ( int k1 = 0; k1 < thrd_num_s; ++k1 )
                        {
                          bucket_count += bucket_count_total[k1];
                          bucket_count_total[k1] = 0;
                        }
                    }

                  if ( io_ready == 2 )
                    {
                      //fprintf(stderr,"All reads have been processed !\n");
                      break;
                    }
                }

              sendWorkSignal ( 3, thrdSignal );
              thread_wait ( threads );
              SwitchBuckets ( &ht2, K_size );

              for ( size_t i = 0; i < ht2.ht_sz; ++i ) //this procedure can be removed
                {
                  struct bucket2 *bktptr = ht2.store_pos[i];

                  while ( bktptr != NULL )
                    {
                      bktptr->kmer_info.cov1 = 0;
                      bktptr = bktptr->nxt_bucket;
                    }
                }

              free ( ( void * ) bucket_count_total );
              free ( ( void * ) locks );
              free ( ( void * ) edge_cnt_total );
              delete [] read_buf0;
              delete [] read_buf1;
              time ( &read_time );
              //fprintf(stderr,"Round 1 consumes time: %.fs.\n",difftime(read_time,beg_time));
              //fprintf(stderr,"Number of processed reads: %llu \n",processed_reads);
              fprintf ( stderr, "Time spent on building graph round 1: %.fs, %llu reads processed, %llu nodes allocated\n",
                        difftime ( read_time, beg_time ), processed_reads, bucket_count );
              fprintf ( stderr, "\n" );
            }

          if ( round == 2 )
            {
              time ( &beg_time );
              edge_cnt = 0;
              //create main io thread
              int read_buf_sz = 102400 * thrd_num_s;
              read_buf0 = new string[read_buf_sz];
              read_buf1 = new string[read_buf_sz];
              io_stat0 = 1; //must be one, if io_stat0 =0 ,the io thread will work immediately
              io_stat1 = 1;
              io_ready = 0;
              io_para_main io_para_mains;
              io_para_mains.read_buf_sz = read_buf_sz;
              io_para_mains.in_filenames_vt = &in_filenames_vt;
              pthread_t io_thread;
              int temp;

              //fprintf(stderr,"Creating main io thread ...\n");
              if ( ( temp = pthread_create ( &io_thread, NULL, run_io_thread_main, &io_para_mains ) ) != 0 )
                {
                  fprintf ( stderr, "ERROR: failed creating main io thread.\n" );
                  exit ( -1 );
                }

              fprintf ( stderr, "1 io thread initialized.\n" );
              //create work threads  for round 2
              pthread_t threads[thrd_num_s];
              unsigned char thrdSignal[thrd_num_s + 1];
              PARAMETER paras[thrd_num_s];
              locks = ( pthread_spinlock_t * ) calloc ( ht2.ht_sz, sizeof ( pthread_spinlock_t ) );

              //initialize locks unlock
              for ( size_t i = 0; i < ht2.ht_sz; ++i )
                {
                  locks[i] = 1;
                }

              //create threads
              //fprintf(stderr,"Creating work threads ...\n");
              bucket_count_total = ( size_t * ) calloc ( thrd_num_s, sizeof ( size_t ) );
              edge_cnt_total = ( size_t * ) calloc ( thrd_num_s, sizeof ( size_t ) );

              for ( int k = 0; k < thrd_num_s; k++ )
                {
                  thrdSignal[k + 1] = 0;
                  paras[k].threadID = k;
                  paras[k].mainSignal = &thrdSignal[0];
                  paras[k].selfSignal = &thrdSignal[k + 1];
                  paras[k].ht = &ht2;
                  paras[k].K_size = K_size;
                  paras[k].gap = gap;
                }

              creatThrds ( threads, paras );
              thrdSignal[0] = 0;
              //begin to work
              int foundcount = 0;
              int flipcount = 0;
              size_t processed_reads = 0;

              while ( 1 )
                {
                  sendIOWorkSignal();

                  while ( io_ready == 0 )
                    {
                      usleep ( 1 );
                    }

                  if ( io_ready )
                    {
                      //read_c = read_num;
                      processed_reads += read_num;
                      sendWorkSignal ( 11, thrdSignal );

                      for ( int k1 = 0; k1 < thrd_num_s; ++k1 )
                        {
                          edge_cnt += edge_cnt_total[k1];
                          edge_cnt_total[k1] = 0;
                        }
                    }

                  if ( io_ready == 2 )
                    {
                      //fprintf(stderr,"All reads have been processed !\n");
                      break;
                    }
                }

              sendWorkSignal ( 3, thrdSignal );
              thread_wait ( threads );
              free ( ( void * ) bucket_count_total );
              free ( ( void * ) locks );
              free ( ( void * ) edge_cnt_total );
              delete [] read_buf0;
              delete [] read_buf1;
              SavingSparseKmerGraph2 ( &ht2, graphfile );
              time ( &read_time );
              //fprintf(stderr,"Round 2 consumed time: %.fs.\n",difftime(read_time,beg_time));
              //fprintf(stderr,"Number of processed reads: %llu \n",processed_reads);
              fprintf ( stderr, "Time spent on building graph round 2: %.fs, %llu reads processed\n", difftime ( read_time, beg_time ), processed_reads );
              fprintf ( stderr, "%llu nodes allocated, %llu kmer-edges allocated.\n", bucket_count, edge_cnt );
              fprintf ( stderr, "\n" );
            }
        }
    }

  if ( LOAD_GRAPH )
    {
      fprintf ( stderr, "Loading the graph ...\n" );
      LoadingSparseKmerGraph2 ( &ht2, graphfile );
    }

  if ( BUILD_EDGES )
    {
      RemovingWeakNodesAndEdges2 ( &ht2, K_size, NodeCovTh, EdgeCovTh, &bucket_count, &edge_cnt );
      int cut_len_tip = 2 * K_size;
      int tip_c = 0;
      time_t start, finish, interval;
      fprintf ( stderr, "Start to remove tips with minority links.\n" );
      start = time ( NULL );
      removeMinorTips ( &ht2, K_size, cut_len_tip, tip_c );
      finish = time ( NULL );
      interval = ( finish - start ) ;
      //fprintf(stderr,"Removing minor tips consumes %llu s.\n\n",interval);
      fprintf ( stderr, "Time spent on removing tips: %llus.\n\n", interval );
      fprintf ( stderr, "Start to construct edges.\n" );
      start = time ( NULL );
      char outfile[256];
      sprintf ( outfile, "%s.sparse.edge", graphfile );
      kmer2edges ( &ht2, K_size, outfile );
      free_hashtable ( &ht2 );
      char temp[256];
      sprintf ( temp, "%s.sparse.edge", graphfile );
      convert ( temp, K_size, graphfile );
      finish = time ( NULL );
      interval = ( finish - start ) ;
      //fprintf(stderr,"Building edges consumes %llu s.\n\n",interval);
      fprintf ( stderr, "Time spent on constructing edges: %llus.\n\n", interval );
    }

  if ( BUILD_PREARCS ) //build preArc
    {
      size_t v_sz, e_sz;
      int K_size;
      char basicInfo[128];
      sprintf ( basicInfo, "%s.preGraphBasic", graphfile );
      FILE *fin;
      char line[1024];
      char str[32];
      fin = fopen ( basicInfo, "r" );

      if ( !fin )
        {
          fprintf ( stderr, "ERROR: can't open file %s\n", basicInfo );
          exit ( 1 );
        }

      bool a = 0, b = 0;

      while ( fgets ( line, 1024, fin ) != NULL )
        {
          if ( line[0] == 'V' ) //VERTEX
            {
              sscanf ( line + 6, "%lu %s %d", &v_sz, str, &K_size );
              a = 1;
            }

          if ( line[0] == 'E' ) //EDGEs
            {
              sscanf ( line, "%s %lu", str, &e_sz );
              b = 1;
              break;
            }
        }

      if ( !a || !b )
        {
          fprintf ( stderr, "ERROR: preGraphBasic file is in invaild format!\n" );
          exit ( 1 );
        }

      vertex_hash2 v_ht;
      preArc_array arc_arr;
      char edge_file[128];
      sprintf ( edge_file, "%s.sparse.edge", graphfile );
      time_t start, finish, interval;
      //step1:
      //fprintf(stderr,"Building vertexes ...\n");
      fprintf ( stderr, "Start to build vertex indexes.\n" );
      start = time ( NULL );
      init_vertex_hash ( &v_ht, v_sz );
      build_vertexes ( &v_ht, K_size, edge_file );
      finish = time ( NULL );
      interval = ( finish - start ) ;
      //fprintf(stderr,"Building vertexes consumes %llu s.\n\n",interval);
      fprintf ( stderr, "Time spent on building vertex indexes: %llus.\n\n", interval );
      //step2:
      //fprintf(stderr,"Building preArcs ...\n");
      fprintf ( stderr, "Start to build preArcs.\n" );
      start = time ( NULL );
      int cut_off_len = 256;//tmp
      init_preArc_array ( &arc_arr, e_sz + 1 );

      if ( solve )
        {
          /* initialize the threads' common data mark_on_edge s_locks*/
          mark_on_edge = ( unsigned int * ) calloc ( e_sz + 1, sizeof ( unsigned int ) );
          s_locks = ( pthread_spinlock_t * ) calloc ( e_sz + 1, sizeof ( pthread_spinlock_t ) );

          for ( size_t i = 0; i < e_sz + 1; i++ )
            {
              s_locks[i] = 1;
            }

          /*buffers for seperate threads*/
          path_buffer = ( edge_path_buffer ** ) calloc ( thrd_num_s, sizeof ( edge_path_buffer ) );

          for ( int i = 0; i < thrd_num_s; i++ )
            {
              path_buffer[i] = create_edge_path_buffer ( mark_on_edge, s_locks, buff_size, max_path_length );
            }

          /*initialize the output file */
          char mark_file[128], path_file[128];
          sprintf ( mark_file, "%s.markOnEdge", graphfile );
          sprintf ( path_file, "%s.path", graphfile );
          mark_fp = fopen ( mark_file, "w" );
          path_fp = fopen ( path_file, "w" );
          pthread_mutex_init ( &file_lock, NULL );
        }

      build_preArc_threaded ( &arc_arr, &v_ht, K_size, cut_off_len, &in_filenames_vt, thrd_num_s );

      if ( solve )
        {
          //output mark_on_edge
          size_t markCounter = 0;

          for ( size_t i = 1; i <= e_sz; i++ )
            {
              markCounter += mark_on_edge[i];
              fprintf ( mark_fp, "%d\n", mark_on_edge[i] );
            }

          fprintf ( stderr, "Total number of marks in file markOnEdge: %lu\n", markCounter );
          fclose ( mark_fp );

          //output path_buffer
          for ( int i = 0; i < thrd_num_s; i++ )
            {
              output_edge_path_buffer ( path_buffer[i], path_fp );
            }

          fclose ( path_fp );

          //destory buffers ...
          for ( int i = 0; i < thrd_num_s; i++ )
            {
              destory_edge_path_buffer ( path_buffer[i] );
            }

          free ( ( void * ) mark_on_edge );
          free ( ( void * ) s_locks );
          pthread_mutex_destroy ( &file_lock );
        }

      finish =  time ( NULL );
      interval = ( finish - start );
      //fprintf(stderr,"Building preArcs consumes %llu s.\n\n",interval);
      fprintf ( stderr, "Time spent on building preArcs: %llus.\n\n", interval );
      char prearc_file[128];
      sprintf ( prearc_file, "%s.preArc", graphfile );
      output_preArcs ( &arc_arr, prearc_file );
    }

  time ( &all_end_time );
  fprintf ( stderr, "Overall time spent on constructing lightgraph: %.fm.\n", difftime ( all_end_time, all_beg_time ) / 60 );
  return 0;
}

static void initenv ( int argc, char **argv )
{
  int copt;
  int inpseq, outseq, genome_sz;
  extern char *optarg;
  char temp[100];
  optind = 1;
  inpseq = outseq = genome_sz = 0;

  while ( ( copt = getopt ( argc, argv, "s:o:K:g:z:d:e:p:m:r:R" ) ) != EOF )
    {
      switch ( copt )
        {
        case 's':
          inpseq = 1;
          sscanf ( optarg, "%s", shortrdsfile );
          continue;

        case 'o':
          outseq = 1;
          sscanf ( optarg, "%s", graphfile );
          continue;

        case 'K':
          sscanf ( optarg, "%s", temp );
          K_size = atoi ( temp );
          continue;

        case 'g':
          sscanf ( optarg, "%s", temp );
          gap = atoi ( temp );
          continue;

        case 'z':
          genome_sz = 1;
          sscanf ( optarg, "%Lu", &GenomeSize );
          continue;

        case 'd':
          sscanf ( optarg, "%s", temp );
          NodeCovTh = atoi ( temp ) >= 0 ? atoi ( temp ) : 0;
          continue;

        case 'e':
          sscanf ( optarg, "%s", temp );
          EdgeCovTh = atoi ( temp ) >= 0 ? atoi ( temp ) : 0;
          continue;

        case 'R':
          solve = 1;
          continue;

        case 'p':
          sscanf ( optarg, "%s", temp );
          thrd_num_s = atoi ( temp );
          continue;

        case 'm':
          continue;

        case 'r':
          sscanf ( optarg, "%s", temp );
          run_mode = atoi ( temp );
          continue;

        default:

          if ( inpseq == 0 || outseq == 0 )
            {
              display_pregraph_usage();
              exit ( -1 );
            }
        }
    }

  if ( inpseq == 0 || outseq == 0 || genome_sz == 0 )
    {
      display_pregraph_usage();
      exit ( -1 );
    }
}

static void parse_args ( vector<string> &in_filenames_vt )
{
  if ( K_size % 2 == 0 )
    {
      K_size--;
    }

#ifdef _63MER_

  if ( K_size > 63 )
    {
      fprintf ( stderr, "ERROR: Parameter K is set too large, max value is 63.\n" );
      exit ( -1 );
    }

#endif
#ifdef _127MER_

  if ( K_size > 127 )
    {
      fprintf ( stderr, "ERROR: Parameter K is set too large, max value is 127.\n" );
      exit ( -1 );
    }

#endif

  if ( gap > 25 )
    {
      fprintf ( stderr, "ERROR: Parameter g is set too large, max value is 25.\n" );
      exit ( -1 );
    }

  //print the args
  fprintf ( stderr, "********************\n" );
  fprintf ( stderr, "SparsePregraph\n" );
  fprintf ( stderr, "********************\n" );
  fprintf ( stderr, "\n" );
  fprintf ( stderr, "Parameters: " );
  fprintf ( stderr, "pregraph_sparse -s %s", shortrdsfile );
  fprintf ( stderr, " -K %d", K_size );
  fprintf ( stderr, " -g %d", gap );
  fprintf ( stderr, " -z %lu", GenomeSize );
  fprintf ( stderr, " -d %d", NodeCovTh );
  fprintf ( stderr, " -e %d", EdgeCovTh );

  if ( solve )
    {
      fprintf ( stderr, " -R " );
    }

  fprintf ( stderr, " -r %d", run_mode );
  fprintf ( stderr, " -p %d", thrd_num_s );
  fprintf ( stderr, " -o %s\n\n", graphfile );

  if ( run_mode == 0 ) //build all
    {
      LOAD_GRAPH = 0;
      BUILD_DBG = 1;
      BUILD_EDGES = 1;
      BUILD_PREARCS = 1;
    }
  else if ( run_mode == 1 ) //build edges ,build preArcs
    {
      LOAD_GRAPH = 1;
      BUILD_EDGES = 1;
      BUILD_PREARCS = 1;
    }
  else if ( run_mode == 2 ) //build graph only
    {
      LOAD_GRAPH = 0;
      BUILD_DBG = 1;
      BUILD_EDGES = 0;
      BUILD_PREARCS = 0;
    }
  else if ( run_mode == 3 ) //build edges only
    {
      LOAD_GRAPH = 1;
      BUILD_EDGES = 1;
      BUILD_PREARCS = 0;
    }
  else if ( run_mode == 4 ) //build preArc only
    {
      LOAD_GRAPH = 0;
      BUILD_DBG = 0;
      BUILD_EDGES = 0;
      BUILD_PREARCS = 1;
    }
  else
    {
      fprintf ( stderr, "ERROR: invalid runMode!\n" );
      exit ( -1 );
    }

  read_lib ( in_filenames_vt, shortrdsfile );
  /*
  fprintf(stderr,"The reads for building sparse de Brujin graph:\n");
  for (int i=0;i<in_filenames_vt.size();++i){
          fprintf(stderr,"%s\n",in_filenames_vt[i].c_str());
  }
  fprintf(stderr,"\n");
  */
}

static void display_pregraph_usage()
{
  fprintf ( stderr, "Usage: sparse_pregraph -s configFile -K kmer -z genomeSize -o outputGraph [-g maxKmerEdgeLength -d kmerFreqCutoff -e kmerEdgeFreqCutoff -R -r runMode -p n_cpu]\n" );
  fprintf ( stderr, "  -s <string>     configFile: the config file of solexa reads\n" );
#ifdef _63MER_
  fprintf ( stderr, "  -K <int>        kmer(min 13, max 63): kmer size, [23]\n" );
#endif
#ifdef _127MER_
  fprintf ( stderr, "  -K <int>        kmer(min 13, max 127): kmer size, [23]\n" );
#endif
  fprintf ( stderr, "  -g <int>        maxKmerEdgeLength(min 1, max 25): number of skipped intermediate kmers, [15]\n" );
  fprintf ( stderr, "  -z <int>        genomeSize(required): estimated genome size\n" );
  fprintf ( stderr, "  -d <int>        kmerFreqCutoff: delete kmers with frequency no larger than,[1]\n" );
  fprintf ( stderr, "  -e <int>        kmerEdgeFreqCutoff: delete kmers' related edge with frequency no larger than [1]\n" );
  fprintf ( stderr, "  -R (optional)   output extra information for resolving repeats in contig step, [NO]\n" );
  fprintf ( stderr, "  -r <int>        runMode: 0 build graph & build edge and preArc, 1 load graph by prefix & build edge and preArc, 2 build graph only, 3 build edges only, 4 build preArcs only [0] \n" );
  fprintf ( stderr, "  -p <int>        n_cpu: number of cpu for use,[8]\n" );
  fprintf ( stderr, "  -o <int>        outputGraph: prefix of output graph file name\n" );
}










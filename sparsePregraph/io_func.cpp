/*
 * io_func.cpp
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
#include "io_func.h"
#include "stdinc.h"

#include <math.h>
#include <unistd.h>

#include "bam.h"
#include "faidx.h"
#include "knetfile.h"
#include "sam_view.h"
#include "xcurses.h"
#include "zlib.h"
#include "bgzf.h"
#include "glf.h"
#include "kstring.h"
#include "razf.h"
#include "sam_header.h"
#include "zconf.h"

#include "sam.h"
#include "errno.h"




//for bam reading ...
static int state = -3;
static int readstate = 0;

static samfile_t *openFile4readb ( const char *fname );

void read1seqbam ( char *src_seq, char *src_name, int max_read_len, samfile_t *in, int *type );


/*************************************************
Function:
    filter_N
Description:
    filters N on reads.
Input:
    1. seq_s:       input read sequence
Output:
    1. seq_s:       read sequence without N
    2. bad_flag:        indicat whether the sequence is vaild
Return:
    None.
*************************************************/
void filter_N ( string &seq_s, int &bad_flag )
{
  //global max_rd_len
  //global min_rd_len
  if ( seq_s.size() > max_rd_len )
    {
      max_rd_len = seq_s.size();
    }

  if ( seq_s.size() < min_rd_len )
    {
      min_rd_len = seq_s.size();
    }

  bad_flag = 0;

  if ( seq_s[seq_s.size() - 1] == '\n' || seq_s[seq_s.size() - 1] == '\r' )
    {
      seq_s.resize ( seq_s.size() - 1 );
    }

  int seq_sz = seq_s.size();
  int nN = seq_sz, isN = -1;

  for ( int i = 0; i < seq_sz; ++i )
    {
      if ( seq_s[i] == '-' || seq_s[i] == 'N' )
        {
          if ( i <= seq_sz / 2 )
            {
              isN = i;
              continue;
            }
          else
            {
              nN = i - 1;
              break;
            }
        }
    }

  if ( nN == seq_sz && isN == -1 )
    {
      bad_flag = 0;
      return;
    }

  if ( ( nN - isN ) <= seq_sz / 2 )
    {
      bad_flag = 1;
    }

  if ( bad_flag == 1 )
    {
      return;
    }

  seq_s = seq_s.substr ( isN + 1, nN - isN );
}



/*************************************************
Function:
    read_lib
Description:
    When the asm_flags equals 1 or 3, adds the filename into vector filenames.
Input:
    1. filenames:       filenames vector
    2. lib_file:        read lib config file
Output:
    1. filenames:       filenames vector
Return:
    None.
*************************************************/
void read_lib ( vector<string> &filenames, char *lib_file )
{
  ifstream lib_in ( lib_file );
  string str;
  int read_stat = 0; // 1: begin a lib, 2:asm_flags=1 or 3
  size_t found;
  int asm_flags;

  while ( getline ( lib_in, str ) )
    {
      if ( read_stat == 0 ) //not start a lib
        {
          found = str.find ( "[LIB]" );

          if ( found == string::npos )
            {
              continue;
            }
          else
            {
              read_stat = 1;
              asm_flags = 0;
            }
        }
      else if ( read_stat == 1 ) //start reading a lib fetch asm flags
        {
          //split by "="
          found = str.find ( "asm_flags" );

          if ( found == string::npos )
            {
              continue;
            }
          else
            {
              found = str.find ( "=" );
              str = str.substr ( found + 1, str.size() - found );

              if ( str.size() == 0 )
                {
                  fprintf ( stderr, "ERROR: please check asm_flags in lib file\n" );
                  exit ( -1 );
                }

              asm_flags = atoi ( str.c_str() );

              if ( asm_flags == 1 || asm_flags == 3 )
                {
                  read_stat = 2;
                }
              else
                {
                  read_stat = 0; // next lib
                }
            }
        }
      else if ( read_stat == 2 ) // reading file
        {
          found = str.find_first_of ( "fqpb" );

          if ( found == 0 ) //f1 f2 q1 q2 p b
            {
              found = str.find ( "=" );

              if ( found > 2 )
                {
                  continue;   // the "=" should be the second or thrid poistion
                }

              str = str.substr ( found + 1, str.size() - found );
              filenames.push_back ( str );
            }
          else
            {
              found = str.find ( "[LIB]" );

              if ( found == string::npos )
                {
                  continue;
                }
              else
                {
                  read_stat = 1;
                  asm_flags = 0;
                }
            }
        }
    }
}


/*************************************************
Function:
    sendIOWorkSignal
Description:
    Initializes the io status, makes the io threads ready to work.
Input:
    None
Output:
    None
Return:
    None
*************************************************/
void sendIOWorkSignal()
{
  if ( io_ready == 2 )
    {
      return ;    //finish io job
    }

  io_stat0 = 0;
  io_stat1 = 0;
  io_ready = 0;
}


/*************************************************
Function:
    run_io_thread_main
Description:
    This is the main io thread working rountine.
    Two buffer "read_buf0" and "read_buf0" are applied to implement AIO.
    It makes the pointer "seq_t" always pointing to the filled buffer.
Input:
    1. arg:     args for io threads
Output:
    None.
Return:
    None.
*************************************************/

void *run_io_thread_main ( void *arg )
{
  io_para_main *paras;
  paras = ( io_para_main * ) arg;

  if ( !paras )
    {
      fprintf ( stderr, "ERROR: the argument passed to main io thread is NULL!\n" );
      exit ( -1 );
    }

  int read_buf_sz = paras->read_buf_sz;
  vector<string> in_filenames_vt = * ( paras->in_filenames_vt );
  /*
  for(int i=0;i<in_filenames_vt.size();i++){
          printf("%s\n",  in_filenames_vt[i].c_str());
  }*/
  int read_num0 = 0;
  int read_num1 = 0;
  char line[1024];
  int read_buf_len = 1024;
  FILE *fp = NULL;  //for normal and gzip reading
  samfile_t *fp3 = NULL;  //for bam reading
  int filetype = 0; //0 normal 1 gzip 2 bam
  int file_num = 0;

  if ( in_filenames_vt.size() >= 1 )
    {
      string temp;
      size_t found;
      found = in_filenames_vt[file_num].find_last_of ( "." );

      if ( found == string::npos )
        {
          //fp = fopen(in_filenames_vt[file_num].c_str(),"r"); //normal
          fp = ( FILE * ) open_file_robust ( "plain", in_filenames_vt[file_num].c_str(), "r" );
          filetype = 0;
        }
      else
        {
          temp = in_filenames_vt[file_num].substr ( found );

          if ( temp.compare ( ".gz" ) == 0 ) //gzip
            {
              //temp = "gzip -dc ";
              //temp.append(in_filenames_vt[file_num]);
              //fp = popen(temp.c_str(),"r");
              fp = ( FILE * ) open_file_robust ( "gz", in_filenames_vt[file_num].c_str(), "r" );
              filetype = 1;
            }
          else if ( temp.compare ( ".bam" ) == 0 ) //bam
            {
              //fp3 = openFile4readb(in_filenames_vt[file_num].c_str());
              fp3 = ( samfile_t * ) open_file_robust ( "bam", in_filenames_vt[file_num].c_str(), "r" );
              filetype = 2;
            }
          else
            {
              //fp = fopen(in_filenames_vt[file_num].c_str(),"r"); //normal
              fp = ( FILE * ) open_file_robust ( "plain", in_filenames_vt[file_num].c_str(), "r" );
              filetype = 0;
            }
        }

      if ( !fp && !fp3 )
        {
          fprintf ( stderr, "ERROR: can't open file %s \n", in_filenames_vt[file_num].c_str() );
          exit ( 1 );
        }
    }
  else
    {
      fprintf ( stderr, "ERROR: input filenames vector is empty! please check the reads config file,option \"asm_flags\" is requried!\n" );
      exit ( 1 );
    }

  //fprintf(stderr,"processing file %d %s \n",file_num,in_filenames_vt[file_num].c_str());
  fprintf ( stderr, "Import reads from file:\n %s\n", in_filenames_vt[file_num].c_str() );

  while ( 1 )
    {
      while ( ( io_stat0 ) && ( io_stat1 ) )
        {
          usleep ( 1 );
        }

      if ( ! ( io_stat0 ) ) //fill buf0
        {
          io_stat0 = 1;// reading reads stat
          int ready = 0;
          int i = 0;

          if ( filetype == 0 || filetype == 1 ) //normal or gzip reading ...
            {
              while ( i < read_buf_sz )
                {
                  if ( fgets ( line, read_buf_len, fp ) != NULL )
                    {
                      switch ( line[0] )
                        {
                        case '@':
                        case '>':
                          ready = 1;
                          break;

                        case '+':
                          ready = 0;
                          break;

                        default:

                          if ( ready )
                            {
                              read_buf0[i].clear();
                              read_buf0[i].append ( line );
                              i++;
                            }
                        }
                    }
                  else
                    {
                      break;
                    }
                }
            }
          else if ( filetype == 2 ) //bam reading
            {
              int type = 0;
              char src_name[128];

              while ( i < read_buf_sz && readstate >= 0 ) //readstate
                {
                  read1seqbam ( line, src_name, read_buf_len, fp3, &type );

                  if ( type != -1 )
                    {
                      read_buf0[i].clear();
                      read_buf0[i].append ( line );
                      //cout<<"line:"<<line<<endl;
                      //cout<<"buf0:"<<read_buf0[i]<<endl;
                      i++;
                    }
                }
            }
          else
            {
              fprintf ( stderr, "ERROR: filetype is not support or filename has a wrong suffix!\n" );
              exit ( -1 );
            }

          read_num0 = i;
          reads_all_num += i;

          while ( io_ready != 0 )
            {
              usleep ( 1 );
            }; //wait for main send work sign

          if ( i == read_buf_sz )
            {
              io_stat0 = 2;
            }
          else if ( i != read_buf_sz && file_num < in_filenames_vt.size() - 1 ) //still has file unread
            {
              io_stat0 = 2;

              if ( filetype == 0 )
                {
                  fclose ( fp );
                }
              else if ( filetype == 1 )
                {
                  pclose ( fp );
                }
              else if ( filetype == 2 )
                {
                  samclose ( fp3 );
                  state = -3;
                  readstate = 0;
                }

              file_num++;
              //open a new file ...
              string temp;
              size_t found;
              found = in_filenames_vt[file_num].find_last_of ( "." );

              if ( found == string::npos )
                {
                  //fp = fopen(in_filenames_vt[file_num].c_str(),"r"); //normal
                  fp = ( FILE * ) open_file_robust ( "plain", in_filenames_vt[file_num].c_str(), "r" );
                  filetype = 0;
                }
              else
                {
                  temp = in_filenames_vt[file_num].substr ( found );

                  if ( temp.compare ( ".gz" ) == 0 ) //gzip
                    {
                      //temp = "gzip -dc ";
                      //temp.append(in_filenames_vt[file_num]);
                      //fp = popen(temp.c_str(),"r");
                      fp = ( FILE * ) open_file_robust ( "gz", in_filenames_vt[file_num].c_str(), "r" );
                      filetype = 1;
                    }
                  else if ( temp.compare ( ".bam" ) == 0 ) //bam
                    {
                      //fp3 = openFile4readb(in_filenames_vt[file_num].c_str());
                      fp3 = ( samfile_t * ) open_file_robust ( "bam", in_filenames_vt[file_num].c_str(), "r" );
                      filetype = 2;
                    }
                  else
                    {
                      //fp = fopen(in_filenames_vt[file_num].c_str(),"r"); //normal
                      fp = ( FILE * ) open_file_robust ( "plain", in_filenames_vt[file_num].c_str(), "r" );
                      filetype = 0;
                    }
                }

              if ( !fp && !fp3 )
                {
                  fprintf ( stderr, "ERROR: can't open file %s \n", in_filenames_vt[file_num].c_str() );
                  exit ( 1 );
                }

              //fprintf(stderr, "processing file %d %s \n",file_num,in_filenames_vt[file_num].c_str());
              fprintf ( stderr, "Import reads from file:\n %s\n", in_filenames_vt[file_num].c_str() );
            }
          else
            {
              io_stat0 = 3;
            }

          seq_t = read_buf0;
          read_num = read_num0;

          if ( io_stat0 == 3 )
            {
              //printf("Io thread's job is finished! all reads: %llu \n",reads_all_num);

              //close the file...
              if ( filetype == 0 )
                {
                  fclose ( fp );
                }
              else if ( filetype == 1 )
                {
                  pclose ( fp );
                }
              else if ( filetype == 2 )
                {
                  samclose ( fp3 );
                  state = -3;
                  readstate = 0;
                }

              io_ready = 2;
              break;
            }

          io_ready = 1;
        }

      if ( ! ( io_stat1 ) ) //fill buf1
        {
          io_stat1 = 1; //reading...
          int ready = 0;
          int i = 0;

          if ( filetype == 0 || filetype == 1 ) //normal or gzip reading ...
            {
              while ( i < read_buf_sz && fp )
                {
                  if ( fgets ( line, read_buf_len, fp ) != NULL )
                    {
                      switch ( line[0] )
                        {
                        case '@':
                        case '>':
                          ready = 1;
                          break;

                        case '+':
                          ready = 0;
                          break;

                        default:

                          if ( ready )
                            {
                              read_buf1[i].clear();
                              read_buf1[i].append ( line );
                              i++;
                            }
                        }
                    }
                  else
                    {
                      break;
                    }
                }
            }
          else if ( filetype == 2 ) //bam reading
            {
              int type = 0;
              char src_name[128];

              while ( i < read_buf_sz && readstate >= 0 ) //readstate
                {
                  read1seqbam ( line, src_name, read_buf_len, fp3, &type );

                  if ( type != -1 )
                    {
                      read_buf1[i].clear();
                      read_buf1[i].append ( line );
                      i++;
                    }
                }
            }
          else
            {
              fprintf ( stderr, "ERROR: filetype is not support or filename has a wrong suffix!\n" );
              exit ( 1 );
            }

          read_num1 = i;
          reads_all_num += i;

          while ( io_ready != 0 )
            {
              usleep ( 1 );
            }; //wait for main send work sign

          if ( i == read_buf_sz && ( fp || fp3 ) )
            {
              io_stat1 = 2;
            }
          else if ( i != read_buf_sz && file_num < in_filenames_vt.size() - 1 ) //still has file unread
            {
              io_stat1 = 2;

              if ( filetype == 0 )
                {
                  fclose ( fp );
                }
              else if ( filetype == 1 )
                {
                  pclose ( fp );
                }
              else if ( filetype == 2 )
                {
                  samclose ( fp3 );
                  state = -3;
                  readstate = 0;
                }

              file_num++;
              //open a new file ...
              string temp;
              size_t found;
              found = in_filenames_vt[file_num].find_last_of ( "." );

              if ( found == string::npos )
                {
                  //fp = fopen(in_filenames_vt[file_num].c_str(),"r"); //normal
                  fp = ( FILE * ) open_file_robust ( "plain", in_filenames_vt[file_num].c_str(), "r" );
                  filetype = 0;
                }
              else
                {
                  temp = in_filenames_vt[file_num].substr ( found );

                  if ( temp.compare ( ".gz" ) == 0 ) //gzip
                    {
                      //temp = "gzip -dc ";
                      //temp.append(in_filenames_vt[file_num]);
                      //fp = popen(temp.c_str(),"r");
                      fp = ( FILE * ) open_file_robust ( "gz", in_filenames_vt[file_num].c_str(), "r" );
                      filetype = 1;
                    }
                  else if ( temp.compare ( ".bam" ) == 0 ) //bam
                    {
                      //fp3 = openFile4readb(in_filenames_vt[file_num].c_str());
                      fp3 = ( samfile_t * ) open_file_robust ( "bam", in_filenames_vt[file_num].c_str(), "r" );
                      filetype = 2;
                    }
                  else
                    {
                      //fp = fopen(in_filenames_vt[file_num].c_str(),"r"); //normal
                      fp = ( FILE * ) open_file_robust ( "plain", in_filenames_vt[file_num].c_str(), "r" );
                      filetype = 0;
                    }
                }

              if ( !fp && !fp3 )
                {
                  fprintf ( stderr, "ERRPR: can't open file %s \n", in_filenames_vt[file_num].c_str() );
                  exit ( 1 );
                }

              //fprintf(stderr,"processing file %d %s \n",file_num,in_filenames_vt[file_num].c_str());
              fprintf ( stderr, "Import reads from file:\n %s\n", in_filenames_vt[file_num].c_str() );
            }
          else
            {
              io_stat1 = 3;
            }

          seq_t = read_buf1;
          read_num = read_num1;

          if ( io_stat1 == 3 )
            {
              //fprintf(stderr,"Io thread's job is finished! all reads: %llu \n",reads_all_num);
              if ( filetype == 0 )
                {
                  fclose ( fp );
                }
              else if ( filetype == 1 )
                {
                  pclose ( fp );
                }
              else if ( filetype == 2 )
                {
                  samclose ( fp3 );
                  state = -3;
                  readstate = 0;
                }

              io_ready = 2;
              break;
            }

          io_ready = 1;
        }
    }

  return NULL;
}


/*************************************************
Function:
    read1seqbam
Description:
    Reads sequence from bam and write it into *src_seq.
Input:
    1. max_read_len:        max read length
    2. in:      sam file
Output:
    1. src_seq:     read sequence
    2. src_name:        read name
        3. type:                record the "state" situation
Return:
    None.
*************************************************/
void read1seqbam ( char *src_seq, char *src_name, int max_read_len, samfile_t *in, int *type )      //read one sequence from bam file
{
  bam1_t *b = bam_init1 ();
  char c;
  char *line1 = NULL;
  int n = 0;
  int len;
  int i, j;
  char *seq1;
  unsigned int flag1 = 0;
  *type = 0;
  readstate = 0;

  if ( ( readstate = samread ( in, b ) ) >= 0 )
    {
      if ( !__g_skip_aln ( in->header, b ) )
        {
          line1 = bam_format1_core ( in->header, b, in->type >> 2 & 3 );
        }

      //printf("%s\n", line2);
      seq1 = strtok ( line1, "\t" );

      for ( i = 0; i < 10; i++ )
        {
          if ( i == 0 )
            {
              sscanf ( seq1, "%s", src_name );
            }
          else if ( i == 1 )
            {
              flag1 = atoi ( seq1 );

              if ( flag1 & 0x0200 )   //whether it's good or not
                {
                  //state(1st read state, 2nd read state) : -3(init), -2(0), -1(1), 0(0, 0), 1(0, 1), 2(1, 0), 3(1, 1)
                  switch ( state )
                    {
                    case -3:
                      state = -2;
                      break;

                    case -2:
                      state = 0;
                      break;

                    case -1:
                      state = 2;
                      break;

                    default:
                      state = -3;
                    }
                }
              else
                {
                  switch ( state )
                    {
                    case -3:
                      state = -1;
                      break;

                    case -2:
                      state = 1;
                      break;

                    case -1:
                      state = 3;
                      break;

                    default:
                      state = -3;
                    }
                }
            }
          else if ( i == 9 )      //the sequence
            {
              //printf("%s\n", seq1);
              len = strlen ( seq1 );

              if ( len + n > max_read_len )
                {
                  len = max_read_len - n;
                }

              for ( j = 0; j < len; j++ )
                {
                  if ( seq1[j] >= 'a' && seq1[j] <= 'z' )
                    {
                      src_seq[n++] = ( char ) ( seq1[j] - 'a' + 'A' );
                    }
                  else if ( seq1[j] >= 'A' && seq1[j] <= 'Z' )
                    {
                      src_seq[n++] = seq1[j];
                      // after pre-process all the symbles would be a,g,c,t,n in lower or upper case.
                    }
                  else if ( seq1[j] == '.' )
                    {
                      src_seq[n++] = 'A';
                    }       // after pre-process all the symbles would be a,g,c,t,n in lower or upper case.
                }

              if ( 3 == state )
                {
                  state = -3;
                }
              else
                {
                  if ( 0 == state || 1 == state || 2 == state )
                    {
                      state = -3;
                      *type = -1;
                    }
                }
            }

          seq1 = strtok ( NULL, "\t" );
        }
    }

  free ( line1 );
  bam_destroy1 ( b );
  src_seq[n++] = '\0';
}

/*************************************************
Function:
    openFile4readb
Description:
    Opens a bam file for reads.
Input:
    1. fname        bam file name
Output:
    None.
Return:
    a samfile pointer
*************************************************/
static samfile_t *openFile4readb ( const char *fname )   //open file to read bam file
{
  samfile_t *in;
  char *fn_list = 0;

  if ( ( in = ( samfile_t * ) samopen ( fname, "rb", fn_list ) ) == 0 )
    {
      fprintf ( stderr, "ERROR: Cannot open %s. Now exit to system...\n", fname );
      return NULL;
      //exit (-1);
    }

  if ( in->header == 0 )
    {
      fprintf ( stderr, "ERROR: Cannot read the header.\n" );
      return NULL;
      //exit (-1);
    }

  return ( in );
}




/*************************************************
Function:
    open_file_robust
Description:
    It opens a file in a "failed-try-again" way.
    It supports plain , gzip and bam format.
Input:
    1. filetype the value enumerate here: "plain", "gz", "bam"
    2. path     the file path
    3. mode     read/write (rw) mode for plain type file, NOTE: gz or bam files are read only
Output:
    None
Return:
    A file pointer with void* type
*************************************************/
void *open_file_robust ( const char *filetype, const char *path, const char *mode )
{
  void *fp = NULL;
  const int max_times = 10;
  const int max_sleep = 60;
  int cur_times = 1;
  int cur_sleep = 1;

  while ( !fp )
    {
      if ( strcmp ( filetype, "plain" ) == 0 )
        {
          if ( access ( path, 0 ) == 0 )
            {
              fp = fopen ( path, mode );
            }
        }
      else if ( strcmp ( filetype, "gz" ) == 0 )
        {
          char tmp[256];
          sprintf ( tmp, "gzip -dc %s", path );

          if ( access ( path, 0 ) == 0 )
            {
              fp = popen ( tmp, "r" );
              /*
              if(fp && feof((FILE*)fp)){ //it's useless for "file not found but popen success"  bug
                  pclose((FILE*)fp);
                  fp = NULL;
              }*/
            }
        }
      else if ( strcmp ( filetype, "bam" ) == 0 )
        {
          if ( access ( path, 0 ) == 0 )
            {
              fp = openFile4readb ( path );
            }
        }

      if ( fp )
        {
          //fprintf(stderr,"%llx \n",fp);
          return fp;
        }
      else
        {
          fprintf ( stderr, "ERROR: open file %s failed!\n", path );
          fprintf ( stderr, "try opening it again after %d seconds\n", cur_sleep );
          sleep ( cur_sleep );
          cur_times ++;
          cur_sleep *= 2;

          if ( cur_sleep >= max_sleep )
            {
              cur_sleep = max_sleep;
            }

          if ( cur_times > max_times )
            {
              fprintf ( stderr, "ERROR: can't open file  %s , now exit system !!!", path );
              exit ( -1 );
              return NULL;
            }
        }
    }
}



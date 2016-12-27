#include "stdinc.h"
#include "newhash.h"
#include "extfunc.h"
#include "global.h"


extern int call_scaffold();
extern int call_align();
extern int call_bundle();
extern int data_prepare();

#define MAPPING 0
#define SCAFF 1
#define BUNDLE 2
#define PREPARE 3
#define POTENT 4
static void usage();
int main(int argc, char **argv)
{
  printf("Mapping & Scaffolding module.\n");

  if(argc == 1)
    {
      usage();
      return 0;
    }

  int c = 0;
  int inpseq, outseq;
  //char optarg[256];
  int mode = -1;

  //char temp[100];
  while((c = getopt(argc, argv, "s:g:p:L:t:i:u:c:P:K:MSBDO")) != EOF)
    {
      switch(c)
        {
        case 'M':
          mode = MAPPING;
          break;

        case 'S':
          mode = SCAFF;
          break;

        case 'B':
          mode = BUNDLE;
          break;

        case 'D':
          mode = PREPARE;
          break;

        case 'O':
          mode = POTENT;
          break;

        case 's':
          inpseq = 1;
          shortrdsfile = (char *)ckalloc(256 * sizeof(char));
          strcpy(shortrdsfile, optarg);
          break;

        case 'g':
          outseq = 1;
          graphfile = (char *)ckalloc(256 * sizeof(char));
          strcpy(graphfile, optarg);
          break;

        case 'p':
          thrd_num = atoi(optarg);
          break;

        case 'L':
          ctg_short = atoi(optarg);
          break;

        case 'P':
          OverlapPercent = atof (optarg);
          break;

        case 't':
          close_threshold = atof (optarg);
          break;

        case 'i':
          ins_size_var = atoi (optarg);
          break;

        case 'u':
          bund_threshold = atoi (optarg);
          break;

        case 'c':
          ctg_file = (char *)ckalloc(256 * sizeof(char));
          strcpy(ctg_file, optarg);
          break;

        case 'K':
          overlaplen = atoi(optarg);
          break;

        case 'h':
          usage();
          break;

        case '?':
          usage();
          exit(1);

        default:
          usage();
          exit(1);
        }
    }

  if(mode == -1)
    {
      usage();
      exit(1);
    }
  else if(mode == MAPPING)
    {
      printf("[%s]Mapping mode selected .\n", __FUNCTION__);

      if(outseq == 0 || inpseq == 0)
        {
          usage();
          exit(1);
        }

      call_align();
    }
  else if(mode == SCAFF)
    {
      printf("[%s]Scaffolding mode selected .\n", __FUNCTION__);

      if(outseq == 0)
        {
          usage();
          exit(1);
        }

      call_scaffold();
    }
  else if(mode == BUNDLE)
    {
      printf("[%s]Bundling mode selected .\n", __FUNCTION__);

      if(outseq == 0)
        {
          usage();
          exit(1);
        }

      call_bundle();
    }
  else if(mode == PREPARE)
    {
      printf("[%s]Data prepare mode selected .\n", __FUNCTION__);

      if(outseq == 0 || ctg_file == NULL)
        {
          usage();
          exit(1);
        }

      data_prepare();
    }
  else if(mode == POTENT)
    {
      printf("[%s]Potential analysis mode selected .\n", __FUNCTION__)	;

      if(outseq == NULL)
        {
          usage();
          exit(1);
        }

      potential();
    }

  return 0;
}

static void usage()
{
  printf("parameters:\n");
  printf("global:\n");
  printf("-s\tLibrary file.\n");
  printf("-g\tPrefix of input files.\n");
  printf("-p\tThreads.\n\n");
  printf("Data prepare mode:\n");
  printf("-D\tEnable this mode.\n");
  printf("-K\tKmer.\n");
  printf("-c\tInput contig file.(can't be name prefix.contig)\n\n");
  printf("Mapping mode:\n");
  printf("-M\tEnable this mode.\n\n");
  printf("Bundling mode.\n");
  printf("-B\tEnable this mode.\n");
  printf("-u\tWeight threshold for outputting bundle file.(default 3)\n\n");
  printf("Potential analysis mode.\n");
  printf("-O\tEnable this mode.\n");
  printf("Scaffolding mode:\n");
  printf("-S\tEnable this mode.\n");
  printf("-L\tthreshold for minimum length of contig(default K+2).\n");
  printf("-P\tOverlap percent threshold for a subgraph(default 0.075).\n");
  printf("-t\tOverlap percent threshold for a PE(default 0.2).\n");
  printf("-i\tOverlap length threshold for remove transitive connect(default 20).\n");
}

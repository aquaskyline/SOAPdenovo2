####################################
#
# Use "make 127mer=1" to make 127mer version
# Use "make 63mer=1" to make 63mer version(default)
#
###################################

CC=             gcc
GCCVERSIONMAJOR := $(shell expr `$(CC) -dumpversion | cut -f1 -d.` \>= 4)
GCCVERSIONMINOR := $(shell expr `$(CC) -dumpversion | cut -f2 -d.` \>= 4)
ifdef debug
CFLAGS=         -O0 -g -fomit-frame-pointer #-static #-mcrc32 -march=core2 -msse4.1 -msse4.2
else
CFLAGS=         -O4 -fomit-frame-pointer #-static #-mcrc32 -march=core2 -msse4.1 -msse4.2
endif
DFLAGS=         
OBJS=		arc.o attachPEinfo.o bubble.o check.o compactEdge.o \
		concatenateEdge.o connect.o contig.o cutTipPreGraph.o cutTip_graph.o \
		darray.o dfib.o dfibHeap.o fib.o fibHeap.o \
		hashFunction.o kmer.o lib.o loadGraph.o loadPath.o \
		loadPreGraph.o localAsm.o main.o map.o mem_manager.o \
		newhash.o node2edge.o orderContig.o output_contig.o output_pregraph.o \
		output_scaffold.o pregraph.o prlHashCtg.o prlHashReads.o prlRead2Ctg.o \
		prlRead2path.o prlReadFillGap.o read2scaf.o readInterval.o stack.o\
		readseq1by1.o scaffold.o searchPath.o seq.o splitReps.o \
		cutTip_graph2.o linearEdge.o kmerhash.o read2edge.o iterate.o
PROG=       SOAPdenovo-63mer
INCLUDES=   -Iinc
SUBDIRS=    .
LIBPATH=    -L/lib64 -L/usr/lib64 -L./inc
LIBS=       -pthread -lz -lm 
EXTRA_FLAGS=

BIT_ERR = 0
ifeq (,$(findstring $(shell uname -m), x86_64 ppc64 ia64))
BIT_ERR = 1
endif

ifdef 127mer
CFLAGS += -DMER127
PROG = SOAPdenovo-127mer
else
CFLAGS += -DMER63
PROG = SOAPdenovo-63mer
endif

ifneq (,$(findstring Linux,$(shell uname)))
EXTRA_FLAGS += -Wl,--hash-style=both
LIBS += -lbam -lrt 
endif

ifneq (,$(findstring Unix,$(shell uname)))
EXTRA_FLAGS += -Wl,--hash-style=both
LIBS += -lbam -lrt
endif

ifneq (,$(findstring Darwin,$(shell uname)))
LIBS += -lbammac
endif

ifneq (,$(findstring $(shell uname -m), x86_64))
CFLAGS += -m64
endif

ifneq (,$(findstring $(shell uname -m), ia64))
CFLAGS += 
endif

ifneq (,$(findstring $(shell uname -m), ppc64))
CFLAGS += -mpowerpc64
endif

.SUFFIXES:.c .o

.c.o:
		@printf "Compiling $<...                             \r"; \
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< || echo "Error in command: $(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $<"

all:		clean $(OBJS)

#SOAPdenovo

.PHONY:all clean install

envTest:
		@test $(BIT_ERR) != 1 || sh -c 'echo "Fatal: 64bit CPU and Operating System required!";false;'
		@test $(GCCVERSIONMAJOR) == 1 || sh -c 'echo "GCC version lower than 4.4.0";false;'
		@test $(GCCVERSIONMINOR) == 1 || sh -c 'echo "GCC version lower than 4.4.0";false;'

SOAPdenovo:	envTest $(OBJS)
		@printf "Linking...                                                \r"
		@$(CC) $(CFLAGS) -o $(PROG) $(OBJS) $(LIBPATH) $(LIBS) $(ENTRAFLAGS)
		@printf "$(PROG) compilation done.\n";

clean:
		@rm -fr gmon.out *.o a.out *.exe *.dSYM $(PROG) *~ *.a *.so.* *.so *.dylib
		@printf "$(PROG) cleaning done.\n";

install:
		@cp $(PROG) ../bin/
		@printf "$(PROG) installed at ../bin/$(PROG)\n"

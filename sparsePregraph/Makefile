CC=            g++ # /opt/blc/gcc-4.5.0/bin/gcc #gcc
ifdef debug
CFLAGS=         -O0 -g -fomit-frame-pointer #-static #-mcrc32 -march=core2 -msse4.1 -msse4.2
else
CFLAGS=         -O4 -fomit-frame-pointer #-static #-mcrc32 -march=core2 -msse4.1 -msse4.2
endif

DFLAGS=         
OBJS=		build_graph.o build_edge.o multi_threads.o \
			build_preArc.o pregraph_sparse.o io_func.o\
			global.o convert_soapdenovo.o
PROG=           
INCLUDES=	-I./inc
SUBDIRS=    . 
LIBPATH=	-L/usr/lib64
LIBS=       -pthread -lz -L./inc
EXTRA_FLAGS=
VERSION =    1.0.3

ifdef 127mer
CFLAGS += -D_127MER_
PROG = pregraph_sparse_127mer.v$(VERSION)

else
CFLAGS += -D_63MER_
PROG = pregraph_sparse_63mer.v$(VERSION)
endif


BIT_ERR = 0
ifeq (,$(findstring $(shell uname -m), x86_64 ppc64 ia64))
BIT_ERR = 1
endif

ifneq (,$(findstring Linux,$(shell uname)))
EXTRA_FLAGS += -Wl,--hash-style=both
LIBS += -lbam
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

.SUFFIXES:.cpp .o

.cpp.o:
		@printf "Compiling $<...                             \r"; \
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< || echo "Error in command: $(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $<"

all:		clean $(OBJS) 
#pregraph_sparse

.PHONY:all clean install

envTest:
		@test $(BIT_ERR) != 1 || sh -c 'echo "Fatal: 64bit CPU and Operating System required!";false;'

pregraph_sparse: clean	envTest $(OBJS)
		@printf "Linking...                                                       \r"
		#@$(CC) $(CFLAGS)$(INCLUDES) -o $(PROG)  $(OBJS) $(LIBPATH) $(LIBS) $(ENTRAFLAGS)
		@printf "$(PROG) compilation done.\n";

clean:
		@rm -fr gmon.out *.o a.out *.exe *.dSYM $(PROG) *~ *.a *.so.* *.so *.dylib
		@printf "$(PROG) cleaning done.\n";

install:
		@cp $(PROG) ../bin/
		@printf "$(PROG) installed at ../bin/$(PROG)\n"

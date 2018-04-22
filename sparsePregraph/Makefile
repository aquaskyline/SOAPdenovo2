CXX ?= g++
ifdef debug
CFLAGS ?= -O0 -g -fomit-frame-pointer -DDEBUG
else
CFLAGS ?= -O3 -fomit-frame-pointer -w
endif

DFLAGS =
OBJS =		build_graph.o build_edge.o multi_threads.o \
			build_preArc.o pregraph_sparse.o io_func.o\
			global.o convert_soapdenovo.o
PROG =
INCLUDES = -I/usr/include/bam -I./inc
SUBDIRS = .
LIBPATH =
LIBS = -pthread -lz
EXTRA_FLAGS=
VERSION = 1.0.3

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

ifneq (,$(findstring $(shell uname -m), x86_64))
CXXFLAGS += -m64
endif

ifneq (,$(findstring $(shell uname -m), ia64))
CXXFLAGS +=
endif

ifneq (,$(findstring $(shell uname -m), ppc64))
CXXFLAGS += -mpowerpc64
endif

.SUFFIXES:.cpp .o

.cpp.o:
		@printf "Compiling $<...                               \r"
		@$(CXX) -c $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $< || echo "Error in command: $(CXX) -c $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $<"

all:		$(OBJS) 
		@printf "$(PROG) objects generated.                    \n"
#pregraph_sparse

.PHONY:all clean

envTest:
		@test $(BIT_ERR) != 1 || sh -c 'echo "Fatal: 64bit CPU and Operating System required!";false;'

pregraph_sparse: clean envTest $(OBJS)
		@printf "Linking...                                     \r"
		@$(CXX) $(CXXFLAGS)$(INCLUDES) -o $(PROG)  $(OBJS) $(LIBPATH) $(LIBS) $(ENTRAFLAGS)
		@printf "$(PROG) compilation done.                      \n"

clean:
		@rm -fr gmon.out *.o a.out $(PROG)
		@printf "$(PROG) cleaning done.                         \n"

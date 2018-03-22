MAKEFLAGS += --no-print-directory
CC = g++
ifdef debug
CFLAGS=         -O0 -g -fomit-frame-pointer -DDEBUG
else
CFLAGS=         -O3 -fomit-frame-pointer
endif

SUBDIRS = sparsePregraph standardPregraph fusion
PROG=       SOAPdenovo-63mer SOAPdenovo-127mer SOAPdenovo-fusion
INCLUDES=   -I./sparsePregraph/inc -I./standardPregraph/inc

LIBPATH=    -L./sparsePregraph/inc -L./standardPregraph/inc -L/lib64 -L/usr/lib64
LIBS=       -pthread -lz -lm
EXTRA_FLAGS= 

BIT_ERR = 0
ifeq (,$(findstring $(shell uname -m), x86_64 ppc64 ia64))
BIT_ERR = 1
endif

ifneq (,$(findstring Linux,$(shell uname)))
EXTRA_FLAGS += -Wl,--hash-style=both
LIBS += -lbam -lrt
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


all: SOAPdenovo-63mer SOAPdenovo-127mer SOAPdenovo-fusion 

SOAPdenovo-fusion:
	@cd fusion;make;cp SOAPdenovo-fusion ../;cd ..;

ifdef debug
SOAPdenovo-63mer:
	@cd sparsePregraph;make 63mer=1 debug=1 clean all;cd ..;
	@cd standardPregraph;make 63mer=1 debug=1 clean all;cd ..;
	@$(CC) sparsePregraph/*.o standardPregraph/*.o $(LIBPATH) $(LIBS) $(EXTRA_FLAGS) -o SOAPdenovo-63mer
SOAPdenovo-127mer:
	@cd sparsePregraph;make 127mer=1 debug=1 clean all;cd ..;
	@cd standardPregraph;make 127mer=1 debug=1 clean all;cd ..;
	@$(CC) sparsePregraph/*.o standardPregraph/*.o $(LIBPATH) $(LIBS) $(EXTRA_FLAGS) -o SOAPdenovo-127mer
else
SOAPdenovo-63mer:
	@cd sparsePregraph;make 63mer=1 clean all;cd ..;
	@cd standardPregraph;make 63mer=1 clean all;cd ..;
	@$(CC) sparsePregraph/*.o standardPregraph/*.o $(LIBPATH) $(LIBS) $(EXTRA_FLAGS) -o SOAPdenovo-63mer
SOAPdenovo-127mer:
	@cd sparsePregraph;make 127mer=1 clean all;cd ..;
	@cd standardPregraph;make 127mer=1 clean all;cd ..;
	@$(CC) sparsePregraph/*.o standardPregraph/*.o $(LIBPATH) $(LIBS) $(EXTRA_FLAGS) -o SOAPdenovo-127mer
endif

clean:
	@cd sparsePregraph;make clean;cd ..;
	@cd standardPregraph;make clean;cd ..;
	@cd fusion;make clean;cd ..;
	@rm -f SOAPdenovo-63mer SOAPdenovo-127mer SOAPdenovo-fusion

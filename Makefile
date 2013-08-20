CC = g++
ifdef debug
CFLAGS=         -O0 -g -fomit-frame-pointer
else
CFLAGS=         -O4 -fomit-frame-pointer
endif

SUBDIRS = sparsePregraph standardPregraph
PROG=       SOAPdenovo-63mer SOAPdenovo-127mer
INCLUDES=   -I./sparsePregraph/inc -I./standardPregraph/inc

LIBPATH=    -L/lib64 -L/usr/lib64 -L./sparsePregraph/inc -L./standardPregraph/inc
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


all: SOAPdenovo-63mer SOAPdenovo-127mer

ifdef debug
SOAPdenovo-63mer:
	@cd sparsePregraph;make 63mer=1 debug=1;cd ..;
	@cd standardPregraph;make 63mer=1 debug=1;cd ..;
	@$(CC) sparsePregraph/*.o standardPregraph/*.o $(LIBPATH) $(LIBS) $(EXTRA_FLAGS) -o SOAPdenovo-63mer
SOAPdenovo-127mer:
	@cd sparsePregraph;make 127mer=1 debug=1;cd ..;
	@cd standardPregraph;make 127mer=1 debug=1;cd ..;
	@$(CC) sparsePregraph/*.o standardPregraph/*.o $(LIBPATH) $(LIBS) $(EXTRA_FLAGS) -o SOAPdenovo-127mer
clean:
	@cd sparsePregraph;make clean;cd ..;
	@cd standardPregraph;make clean;cd ..;
	@rm SOAPdenovo-63mer SOAPdenovo-127mer -f
else
SOAPdenovo-63mer:
	@cd sparsePregraph;make 63mer=1;cd ..;
	@cd standardPregraph;make 63mer=1;cd ..;
	@$(CC) sparsePregraph/*.o standardPregraph/*.o $(LIBPATH) $(LIBS) $(EXTRA_FLAGS) -o SOAPdenovo-63mer
SOAPdenovo-127mer:
	@cd sparsePregraph;make 127mer=1;cd ..;
	@cd standardPregraph;make 127mer=1;cd ..;
	@$(CC) sparsePregraph/*.o standardPregraph/*.o $(LIBPATH) $(LIBS) $(EXTRA_FLAGS) -o SOAPdenovo-127mer
clean:
	@cd sparsePregraph;make clean;cd ..;
	@cd standardPregraph;make clean;cd ..;
	@rm SOAPdenovo-63mer SOAPdenovo-127mer -f
endif

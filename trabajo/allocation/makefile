#####################################################
# Modified 2017-04-03 --jha & xh
# Last Modified 2018-08-30 --ltz1
####################################################
# SRC File name
SRCFILE = solvers_w

HEADERS = 

INCLUDE =

LIBS =

# Switch between debug or optimized
CFLAGS = -g
FFLAGS = -g

WITH_SUITESPARSE=1

# Extension for the Executable Programs
EXTENSION = ex

# Machine Specific Compilers and Libraries
CC = gcc
FC = gfortran
CXX = g++
CFLAGS += 
FFLAGS += -fno-second-underscore

ExtraFLAGS =
INCLUDE =
LIBS =


##################### no change should be needed below. ###########
# HAZMATH LIB and INCLUDE
HAZDIR = $(realpath ../..)

HAZLIB = -L$(HAZDIR)/lib -lhazmath

INCLUDE += -I$(HAZDIR)/include

LIBS += $(HAZLIB) -lm -lblas -llapack

HEADERS += 

INCLUDESSP=
ifeq ($(WITH_SUITESPARSE),1)
	CFLAGS += -DWITH_SUITESPARSE=1
	SSDIR = /usr/lib/x86_64-linux-gnu
	INCLUDESSP = -I/usr/include/suitesparse
	LIBS += -lsuitesparseconfig -lcholmod -lamd -lcolamd -lccolamd -lcamd -lspqr -lumfpack -lamd -lcxsparse 
endif


############### 
# Different Executable Programs, but same targets; SRC file needs to be defined

EXE = $(SRCFILE).$(EXTENSION)

# Source and Object Files
OBJS += $(SRCFILE).o

HEADERS += $(HEADERS)

LIBS += #-lgfortran

.PHONY:  all

all: $(EXE) 

$(EXE):	$(MGTARGET)	$(OBJS)	
	+$(CC) $(CFLAGS) $(ExtraFLAGS) $(INCLUDE) $(OBJS) $(MGLIBS) -o $@  $(LIBS)

%.o:	%.c
	+$(CC) $(INCLUDE) $(INCLUDESSP) $(CFLAGS) $(DMGRAPH) -o $@ -c $<

clean:
	+rm -rf $(EXE) $(OBJS) *.mod output/* ./*.dSYM  $(EXTRA_DEL)




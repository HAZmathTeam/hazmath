#####################################################
# Last Modified 2017-03-08 --ltz
####################################################

# Extension for the Executable Programs
EXTENSION = ex
# Machine Specific Compilers and Libraries
CC = gcc-6
FC = gfortran-6
CXX = g++-6
CFLAGS = -g 
FFLAGS = -g -fno-second-underscore
ExtraFLAGS =
# SuiteSparse  directory with library
#JA (Macports)
#SSDIR = /opt/local
# Polaris 
#SSDIR = /usr/local/numerics/suitesparse
#XH & LZ (Linux | Homebrew)

SSDIR = /usr/local

##################### no change should be needed below. ###########

# HAZMATH LIB and INCLUDE
HAZDIR = $(realpath ../..)
HAZLIB = -L$(HAZDIR)/lib -lhazmath
SSLIB = -L$(SSDIR)/lib -lsuitesparseconfig -lcholmod -lamd -lcolamd -lccolamd -lcamd -lspqr -lumfpack -lamd -lcxsparse 

LIBS = $(HAZLIB) -lm -lblas -llapack -lgfortran $(SSLIB) 

INCLUDE = -I$(HAZDIR)/include  -I$(SSDIR)/include

############### 
# Different Executable Programs, but same targets; SRC file needs to be defined

EXE = $(SRCFILE).$(EXTENSION)

# Source and Object Files
OBJS = $(OBJS_ADD) $(SRCFILE).o
HEADERS=$(HEADERS_ADD)

.PHONY: all

all: $(EXE)

$(EXE):	$(OBJS) $(HEADERS)
	+$(CC) $(ExtraFLAGS) $(INCLUDE) $(OBJS) -o $(EXE)  $(LIBS)

%.o:	%.c $(HEADERS)
	+$(CC) $(INCLUDE) $(CFLAGS) -o $@ -c  $< 

clean:
	+rm -rf $(EXE) $(OBJS) output/* ./*.dSYM  $(EXTRA_DEL)

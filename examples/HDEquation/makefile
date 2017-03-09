#####################################################
# Last Modified 2017-02-11 --jha
####################################################

# SRC File name
SRCFILE = HDEquation

# Different Executable Programs
EXE = $(SRCFILE).ex

# Machine Specific Compilers and Libraries
CC = gcc
FC = gfortran
CXX = g++
CFLAGS = -g 
FFLAGS = -g -fno-second-underscore
ExtraFLAGS =

# HAZMATH
HAZDIR=../..
HAZLIB = $(HAZDIR)/lib
HAZINCLUDE = -I$(HAZDIR)/include

# SuiteSparse  
#XH (Linux)
#SSDIR = /usr/local
#JA (Macports)
SSDIR = /opt/local
#LTZ (Homebrew)
SSDIR = /usr/local
# Polaris 
#SSDIR = /usr/local/numerics/suitesparse
SSINCLUDE = -I$(SSDIR)/include
SSLIB = -L$(SSDIR)/lib -lsuitesparseconfig -lcholmod -lamd -lcolamd -lccolamd -lcamd -lspqr -lumfpack -lamd -lcxsparse 

# Source and Object Files
SRC = .
SRCex = $(SRC)/$(SRCFILE).c
OBJex = $(SRC)/$(SRCFILE).o
OBJSex = $(OBJex) $(HAZLIB)/libhazmath.a 

# Include all Libraries including UMFPACK, its dependencies, BLAS, and LAPACK
INCLUDE = $(SSINCLUDE) $(HAZINCLUDE)
LIBS = -lm -lblas -llapack $(SSLIB)

.PHONY: all

all: $(EXE)

$(EXE): $(OBJSex)
	$(CC) $(ExtraFLAGS) $(OBJSex) $(INCLUDE) $(LIBS) -o $(EXE) 

%.o:	%.c
	$(CC) $(INCLUDE) $(CFLAGS) -o $@ -c  $< $(PROFFLAG)

%.o:	%.f
	$(FC) $(INCLUDE) $(FFLAGS) -o $@ -c  $< $(PROFFLAG)

clean:
	rm -rf $(EXE) $(SRC)/*.o output/* ./*.dSYM


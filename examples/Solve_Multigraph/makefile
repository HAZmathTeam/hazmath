#####################################################
# Last Modified 2017-03-06 --ltz
####################################################

# SRC File name
SRCFILE = mgraph

# Different Executable Programs
EXE = $(SRCFILE).ex

# Machine Specific Compilers and Libraries
CC = gcc

# HAZMATH
HAZDIR=../..
HAZLIB = $(HAZDIR)/lib
HAZINCLUDE = -I$(HAZDIR)/include

# SuiteSparse  
# XH (Linux)
#SSDIR = /usr/local
#JA (Macports)
SSDIR = /opt/local
# Polaris 
#SSDIR = /usr/local/numerics/suitesparse
SSINCLUDE = -I$(SSDIR)/include
SSLIB = -L$(SSDIR)/lib -lsuitesparseconfig -lcholmod -lamd -lcolamd -lccolamd -lcamd -lspqr -lumfpack -lamd -lcxsparse 

# Source and Object Files
SRC = .
SRCex = $(SRC)/$(SRCFILE).c
OBJex = $(SRC)/$(SRCFILE).o
OBJSex = $(OBJex) $(HAZLIB)/libhazmat.a  $(HAZLIB)/libmultigraph.a 

# Include all Libraries including UMFPACK, its dependencies, BLAS, and LAPACK
INCLUDE = $(SSINCLUDE) $(HAZINCLUDE)
LIBS = -lm -lblas -llapack -lgfortran $(SSLIB)

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


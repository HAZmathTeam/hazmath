#####################################################
# Last Modified 2012-07-21 --jha
####################################################

# Different Executable Programs
EXE = HeatEqn.ex

# Generic Directories
DIR0 = .
#XH
#SSDIR = /usr/local
#JA
#SSDIR = /opt/local
#Polaris
SSDIR = /usr/local/numerics/suitesparse

# Machine Specific Compilers and Libraries
# APPLE MACHINES
CC = gcc
FC = gfortran
CXX = g++
CFLAGS = -g 
FFLAGS = -g -fno-second-underscore
LDFLAGS = -g
SSINCLUDE = -I$(SSDIR)/include
SSLIB = -L$(SSDIR)/lib -lsuitesparseconfig -lcholmod -lamd -lcolamd -lccolamd -lcamd -lspqr -lumfpack -lamd -lcxsparse 

# HAZMAT Location
HAZDIR=../..
LIB = $(HAZDIR)/lib
HAZINCLUDE = -I$(HAZDIR)/include

INCLUDE = $(SSINCLUDE) $(HAZINCLUDE)

SRC = .
SRCex = $(SRC)/HeatEqn.c

OBJex = $(SRC)/HeatEqn.o

OBJS= $(LIB)/libhazmat.a	

OBJSex = $(OBJex) $(OBJS) 

# Include all Libraries including UMFPACK, its dependencies, BLAS, and LAPACK
LIBS = -lm -lblas -llapack -lgfortran $(SSLIB)

.PHONY: all

all: $(EXE)

$(EXE): $(OBJSex)
	$(CC) $(LDFLAGS) $(OBJSex) $(INCLUDE) $(LIBS) -o $(EXE) 

%.o:	%.c
	$(CC) $(INCLUDE) $(CFLAGS) -o $@ -c  $< $(PROFFLAG)

%.o:	%.f
	$(FC) $(INCLUDE) $(FFLAGS) -o $@ -c  $< $(PROFFLAG)

clean:
	rm -rf $(EXE) $(SRC)/*.o sol.dat ./*.dSYM


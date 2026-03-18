#####################################################
# Last Modified 2016-08-30 --jha
####################################################

# Different Executable Programs
EXE = LC.exe
SRCFILE = lcelastic

# Who is compiling?
USER := $(shell logname)

# Machine Specific Compilers and Libraries
CC = gcc
FC = gfortran
CXX = g++
CFLAGS = -g
FFLAGS = -g -fno-second-underscore
LDFLAGS = -g 

# JHA (default)
# HAZ Directory
# Hazmat Directory
#JA
HAZDIR = /Users/jadler05/Documents/James/Research/CODE/HAZMATH/main
#XH
#HAZDIR = /home/xiaozhehu/Work/Projects/HAZMATH/hazmath
#DE
#HAZDIR = /Users/David/Desktop/ResearchRepos/HAZMath/hazmath

#JA
SSDIR = /opt/homebrew
#DE
#SSDIR = /usr/local/
SSDIRlib = $(SSDIR)/lib

# LTZ
ifeq ($(USER),ltz1)
  # HAZ Directory
  HAZDIR = /home/$(USER)/hazmath-stable
  # Suitesparse
  SSDIR = /usr
  SSDIRlib = /usr/lib/x86_64-linux-gnu
endif


# Include all other Libraries including Suitesparse, its dependencies, BLAS, and LAPACK
SSINCLUDE = -I$(SSDIR)/include
SSLIB = -L$(SSDIRlib) -lsuitesparseconfig -lcholmod -lamd -lcolamd -lccolamd -lcamd -lspqr -lumfpack -lcxsparse 
LIBS = -lm -lhazmath -lblas -llapack $(SSLIB) 

# HAZMATH Files
LIB = $(HAZDIR)/lib
HAZINCLUDE = -I$(HAZDIR)/include
INCLUDE = $(SSINCLUDE) $(HAZINCLUDE)

# System specific include and libraries
SYSNAME := $(shell uname -s)
# MAC?
ifeq ($(SYSNAME),Darwin)
	MAC_INCLUDE=-I/opt/homebrew/include
	MAC_LIB=-L/opt/homebrew/lib
	INCLUDE += $(MAC_INCLUDE)
	LIBS += $(MAC_LIB)
	MAC_ON=1
	RPATH = 
else
	MAC_ON=0
	RPATH = -Wl,-rpath=$(HAZDIR)/lib
endif

# Source and Object Files
SRCDIR = .
OBJ = $(SRCDIR)/$(SRCFILE).o
OBJS= $(OBJ) 

.PHONY: all

all: $(EXE)

$(EXE): $(OBJS)
	$(CC)  $(OBJS) -L$(LIB) $(RPATH) $(LIBS) -o $(EXE) 

%.o:	%.c
	$(CC) $(INCLUDE) $(CFLAGS) -o $@ -c  $< $(PROFFLAG)

%.o:	%.f
	$(FC) $(INCLUDE) $(FFLAGS) -o $@ -c  $< $(PROFFLAG)

clean:
	rm -rf $(EXE) $(SRCDIR)/*.o *.dSYM

veryclean:
	rm -rf $(EXE) $(SRCDIR)/*.o output/* *.dSYM


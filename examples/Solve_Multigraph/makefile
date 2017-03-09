#####################################################
# Last Modified 2017-03-06 --ltz
####################################################

SRCDIR = .
# SRC File name
SRCFILE = $(SRCDIR)/multigraph_solve

CC = gcc-6
FC = gfortran-6
CFLAGS = -g
# Different Executable Programs
EXE = $(SRCFILE).ex

# Machine Specific Compilers and Libraries
# HAZMATH
HAZDIR=$(realpath ../..)
HAZLIB = $(HAZDIR)/lib
HAZINCLUDE = $(HAZDIR)/include
#Multigraph
MULTIGRAPH_DIR = $(realpath /Users/ltz1/packages/multigraph_2.0)
##########################################################
MULTIGRAPH_SRCDIR = $(MULTIGRAPH_DIR)/src
MULTIGRAPH_PRECISION_MOD = $(MULTIGRAPH_SRCDIR)/mthdef.mod
MULTIGRAPH_PRECISION_SRC = $(MULTIGRAPH_SRCDIR)/mg0.f


# Source and Object Files
SRCex = $(SRCFILE).c
HEADERS = $(SRCFILE).h

OBJS = solver.o $(SRCFILE).o m000.o
INCLUDE = -I$(HAZINCLUDE) -I$(SRCDIR)
LIBS = -lhazmath -lgfortran


.PHONY: all clean

all: $(EXE)

$(EXE):	$(OBJS)  $(HAZLIB)
	+$(CC)  -L$(HAZLIB) $(OBJS) -o $@ $(LIBS)

%.o:	 %.c $(HEADERS)
	+$(CC)  $(CFLAGS) $(INCLUDE)  -L$(HAZLIB) -o $@ -c $< 

$(MULTIGRAPH_PRECISION_MOD):	$(MULTIGRAPH_PRECISION_SRC)
		+$(FC)  $(CFLAGS) -c $<  $(INCLUDE) -o $@ 

%.o:	 $(MULTIGRAPH_SRCDIR)/%.f $(MULTIGRAPH_PRECISION_MOD)
	+$(FC)  $(CFLAGS) -c $<  $(INCLUDE) -o $@ 

# Generate an error message if the HAZmath library is not built
$(HAZLIB):
	$(error The HAZmath library is not built or is not up to date)

clean:
	rm -f *.o *.mod *~ $(MULTIGRAPH_PRECISION_MOD)  $(EXE) 

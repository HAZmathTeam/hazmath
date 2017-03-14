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

##################### no change should be needed below. ###########
# HAZMATH LIB and INCLUDE
HAZDIR = $(realpath ../..)

HAZLIB = -L$(HAZDIR)/lib -lhazmath

LIBS = $(HAZLIB) -lm -lblas -llapack -lgfortran $(LIBS_ADD)

INCLUDE = -I$(HAZDIR)/include  -I$(SSDIR)/include

ifeq ($(WITH_MGRAPH),yes)
	DMGRAPH = -DMGRAPH
	MGTARGET=multigraph
	MGRAPH_SRCDIR = $(realpath ../../../../multigraph_2.0/src)
	MGRAPH_WRAPPERDIR = $(realpath ../with_multigraph)
	MGLIBS = $(MGRAPH_WRAPPERDIR)/multigraph_solve.o $(MGRAPH_SRCDIR)/solver.o
else
	DMGRAPH =
	MGTARGET = $(EXE)
	MGRAPH_SRCDIR = 
	MGRAPH_WRAPPERDIR = 
	MGLIBS = 
endif

############### 
# Different Executable Programs, but same targets; SRC file needs to be defined

EXE = $(SRCFILE).$(EXTENSION)

# Source and Object Files
OBJS = $(SRCFILE).o

HEADERS = $(HEADERS_ADD)

.PHONY:  all

all: $(EXE) 

$(EXE):	$(MGTARGET)	$(OBJS)
	+$(CC) $(ExtraFLAGS) $(INCLUDE) $(DMGRAPH) $(OBJS) $(MGLIBS) -o $(EXE)  $(LIBS)

%.o:	%.c
	+$(CC) $(INCLUDE) $(CFLAGS) -o $@ -c $<

clean:
	+rm -rf $(EXE) $(OBJS) output/* ./*.dSYM  $(EXTRA_DEL)

multigraph:	
	@if [  -f  $(MGRAPH_SRCDIR)/mg0.f \
		-a  -f $(MGRAPH_SRCDIR)/solver.f \
		-a  -f $(MGRAPH_WRAPPERDIR)/multigraph_solve.c \
		-a  -f $(MGRAPH_WRAPPERDIR)/multigraph_solve.h ] ; then \
		 make -B CC=$(CC) FC=$(FC) $(INCLUDE) \
		$(MGRAPH_SRCDIR)/mg0.o	\
		$(MGRAPH_SRCDIR)/solver.o \
		$(MGRAPH_WRAPPERDIR)/multigraph_solve.o ; \
	else \
		echo "*******************************************************************"; \
		echo "* One or more of the MULTIGRAPH sources were not found. " \
		echo "  MULTIGRAPH directory: $(MULTIGRAPH_SRCDIR) "; \
		echo " Setting WITH_MULTIGRAPH=0"; \
		echo "*******************************************************************"; \
	fi


# Generate an error message if the HAZmath library is not built
# $(HAZLIBFILE):
#	$(error The HAZmath library is not built or is not up to date)


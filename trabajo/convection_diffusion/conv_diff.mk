#####################################################
# Last Modified 2017-03-08 --ltz
####################################################

# Extension for the Executable Programs
EXTENSION = ex
# Machine Specific Compilers and Libraries
CC = gcc
FC = gfortran
CXX = g++
CFLAGS += 
FFLAGS += -fno-second-underscore
ExtraFLAGS =

##################### no change should be needed below. ###########
# HAZMATH LIB and INCLUDE
HAZDIR = $(realpath ../..)

HAZLIB = -L$(HAZDIR)/lib -lhazmath

INCLUDE = -I$(HAZDIR)/include $(INCLUDEADD)

LIBS = $(HAZLIB) -lm -lblas -llapack $(LIBSADD) #-lgfortran

HEADERS += 

MGRAPH_WRAPPERDIR = $(realpath ../../examples/multigraph_wrap)
DMGRAPH = 
ifeq ($(WITH_MGRAPH),yes)
ifneq "$(wildcard $(MGRAPH_SRCDIR) )" ""
ifneq "$(wildcard $(MGRAPH_WRAPPERDIR) )" ""
	DMGRAPH := -DMGRAPH
else
$(warning *** No Multigraph support MGRAPH_WRAPPERDIR="$(MGRAPH_WRAPPERDIR)" ***)
endif
else
$(warning *** No Multigraph support MGRAPH_SRCDIR="$(MGRAPH_SRCDIR)" ***)
endif
endif

ifeq ($(DMGRAPH),) 
	DMGRAPH := -UMGRAPH
	MGTARGET = 
	MGRAPH_SRCDIR = 
	MGRAPH_WRAPPERDIR = 
	MGLIBS = 
else
	MGTARGET := multigraph
	MGLIBS = $(MGRAPH_WRAPPERDIR)/multigraph_solve.o $(MGRAPH_SRCDIR)/solver.o  -lgfortran
endif

############### 
# Different Executable Programs, but same targets; SRC file needs to be defined

EXE = $(SRCFILE).$(EXTENSION)

# Source and Object Files
OBJS += $(SRCFILE).o

HEADERS += $(HEADERS)

.PHONY:  all

all: $(EXE) 

$(EXE):	$(MGTARGET)	$(OBJS)	
	+$(CC) $(ExtraFLAGS) $(INCLUDE) $(OBJS) $(MGLIBS) -o $@  $(LIBS)

%.o:	%.c
	+$(CC) $(INCLUDE) $(CFLAGS) $(DMGRAPH) -o $@ -c $<

%.o:	%.f
	+$(FC) $(INCLUDE) $(FFLAGS) $(DMGRAPH) -o $@ -c $<


clean:
	+rm -rf $(EXE) $(OBJS) *.mod output/* ./*.dSYM  $(EXTRA_DEL)

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
		echo "  MULTIGRAPH directory: $(MGRAPH_SRCDIR) "; \
		echo " Setting MGRAPH=0"; \
		echo "*******************************************************************"; \
	fi


# Generate an error message if the HAZmath library is not built
# $(HAZLIBFILE):
#	$(error The HAZmath library is not built or is not up to date)

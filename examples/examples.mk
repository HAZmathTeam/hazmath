
####################################################
# Modified 2017-03-08 --ltz
# Last Modified 2023-02-02 --ltz
####################################################

# Extension for the Executable Programs
ifeq ($(MAKE_LIB),1)
## making shared lib from a source. 
	CFLAGS += -fPIC -shared
	EXTENSION = so
else
	EXTENSION = ex
endif

# Machine Specific Compilers and Libraries
CC = gcc
FC = gfortran
CXX = g++
CFLAGS += 
FFLAGS += -fno-second-underscore
ExtraFLAGS =
INCLUDE = 
LIBS =

SYSNAME := $(shell uname -s)
# MAC_ON=0 #[0/1] this should be given in the particular makefile.
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

# PROCNAME := $(shell uname -p)
# ifeq ($(PROCNAME),x86_64)
#     CCFLAGS += -D AMD64
# endif
# ifneq ($(filter %86,$(PROCNAME)),)
#     CCFLAGS += -D IA32
# endif
# ifneq ($(filter arm%,$(PROCNAME)),)
#     CCFLAGS += -D ARM
# endif
# MAC specific (assuming homebrew)

##################### no change should be needed below. ###########
# HAZMATH LIB and INCLUDE
ifndef COMPILER
## uncomment this below to see the annoying warnings
#$(warning COMPILER not set; setting COMPILER to $(CC))
COMPILER = $(CC)
endif
ifndef LINKER
## uncomment this below to see the annoying warnings
#$(warning LINKER not set; setting LINKER to $(COMPILER))
LINKER = $(COMPILER)
endif

HAZDIR = $(realpath ../..)

HAZLIB = -L$(HAZDIR)/lib -lhazmath

INCLUDE += -I$(HAZDIR)/include

LIBS += $(HAZLIB) -lm 

MGRAPH_WRAPPERDIR = $(realpath ../multigraph_wrap)
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
	MGLIBS = $(MGRAPH_WRAPPERDIR)/multigraph_solve.o $(MGRAPH_SRCDIR)/solver.o
endif

ifeq ($(WITH_BLAS),1)
	CFLAGS += -DWITH_BLAS=1
	LIBS += -lblas 
endif

ifeq ($(WITH_LAPACK),1)
	CFLAGS += -DWITH_LAPACK=1
	LIBS += -llapack
endif

ifeq ($(WITH_HDF5),1)
	CFLAGS += -DWITH_HDF5=1
	LIBS += -lhdf5_serial
	INCLUDE += -I/usr/include/hdf5/serial
endif

ifeq ($(WITH_FFTW3),1)
	CFLAGS += -DWITH_FFTW3=1
	LIBS += -lfftw3
endif

ifeq ($(WITH_HAZNICS),1)
	CFLAGS += -DWITH_HAZNICS=1
	INCLUDE += -I/usr/include/python3.8
endif

INCLUDESSP=
ifeq ($(WITH_SUITESPARSE),1)
	CFLAGS += -DWITH_SUITESPARSE=1
	ifeq ($(MAC_ON),0)
		SSDIR = /usr/lib/x86_64-linux-gnu
		INCLUDESSP = -I/usr/include/suitesparse
	endif
	LIBS += -lsuitesparseconfig -lcholmod -lamd -lcolamd -lccolamd -lcamd -lspqr -lumfpack -lamd -lcxsparse 
endif


############### 
# Different Executable Programs, but same targets; SRC file needs to be defined

EXE = $(SRCFILE).$(EXTENSION)

# Source and Object Files
OBJS += $(SRCFILE).o

HEADERS += 

LIBS += #-lgfortran

.PHONY:  all

all: $(EXE) 

$(EXE):	$(MGTARGET)	$(OBJS)	
	+$(LINKER) $(CFLAGS) $(ExtraFLAGS) $(OBJS) $(RPATH) $(MGLIBS) -o $@  $(LIBS)

%.o:	%.c	$(HEADERS)
	+$(COMPILER) $(INCLUDE) $(INCLUDESSP) $(CFLAGS) $(DMGRAPH) -o $@ -c $<

%.o:	%.cpp	$(HEADERS)
	+$(COMPILER) $(INCLUDE) $(INCLUDESSP) $(CFLAGS) $(DMGRAPH) -o $@ -c $<

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


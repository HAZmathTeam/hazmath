#####################################################
# Modified 2012-07-21 --jha
# Modified 2015-01-10 --Xiaozhe 
####################################################

########################################################################
# Generic Directories
########################################################################
DIR0 = .
INCLUDE = -I$(DIR0)/include
CSRCDIR = $(DIR0)/src
LIB = $(DIR0)/lib/libHAZMAT.a
DRIVERS = $(DIR0)/drivers

########################################################################
# Machine Specific Compilers
########################################################################
CC = gcc
FC = gfortran
CXX = g++
CPP = g++
AR = ar ruc

########################################################################
# Depedencs 
########################################################################
# ARPACK
#ARPACKLIB = 

# UMFPACK
#XH 
#UMFPACKDIR = /usr/local
#JA
UMFPACKDIR = /opt/local
UMFPACKINCLUDE = -I$(UMFPACKDIR)/include
UMFPACKLIB = -L$(UMFPACKDIR)/lib -lsuitesparseconfig -lcholmod -lamd -lcolamd -lccolamd -lcamd -lspqr -lumfpack -lamd -lcxsparse

# METIS
#METISDIR = /MFEM-SVN/metis-4.0
#METISINCLUDE = -I$(METISDIR)/include
#METISLIB = -L$(METISDIR) -lmetis

# BLAS
BLASLIB = -framework Accelerate

########################################################################
# Compiling Options
########################################################################
BOPT = -g -pg -O0 -Wall #-fopenmp

COPTS = $(BOPT)
CINCLUDES = $(INCLUDE) $(UMFPACKINCLUDE) $(METISINCLUDE)
CFLAGS = $(COPTS) $(CINCLUDES)

FOPTS = $(BOPT) -fno-second-underscore
FINCLUDES = $(CINCLUDES)
FFLAGS = $(FOPTS) $(FINCLUDES)

########################################################################
# Set Libraries 
########################################################################
# Include all Libraries including UMFPACK, its dependencies, BLAS, and LAPACK
LIBS = $(LIB) $(BLASLIB) $(UMFPACKLIB) $(ARPACKLIB) $(METISLIB) 

########################################################################
# Link options
########################################################################
LINKOPTS = $(BOPT)
CLFLAGS = -lstdc++ $(LINKOPTS) $(LIBS)
FLFLAGS = -lm $(LINKOPTS) $(LIBS)

########################################################################
# Rules for compiling the source files
########################################################################
.SUFFIXES: .c .cc .cpp
#
CSRC := $(wildcard $(CSRCDIR)/assemble/*.c) 
CSRC += $(wildcard $(CSRCDIR)/fem/*.c)
CSRC += $(wildcard $(CSRCDIR)/grid/*.c)
CSRC += $(wildcard $(CSRCDIR)/solver/*.c)
CSRC += $(wildcard $(CSRCDIR)/utilities/*.c)
#
OBJSC := $(patsubst %.c,%.o,$(CSRC))
#
.c.o:
	$(CC) -c $< -o $@ $(CFLAGS) 
	$(AR) $(LIB) $@
#
.cc.o:
	$(CPP) -c $< -o $@ $(CFLAGS) 
	$(AR) $(LIB) $@
#
.cpp.o:
	$(CPP) -c $< -o $@ $(CFLAGS) 
	$(AR) $(LIB) $@
#

########################################################################
# List of all programs to be compiled
########################################################################

# Everything
ALLPROG=$(LIB) test diff

########################################################################
# Link
########################################################################

all: $(ALLPROG)

Default: $(ALLPROG)	

headers: 
	/bin/cat $(CSRCDIR)/assemble/*.c $(CSRCDIR)/fem/*.c \
	 $(CSRCDIR)/grid/*.c $(CSRCDIR)/solver/*.c  $(CSRCDIR)/utilities/*.c \
	| awk -v name="functs.h" -f mkheaders.awk > ./include/functs.h

$(LIB): $(OBJSC)
	ranlib $(LIB)

lib: $(OBJSC)
	ranlib $(LIB)

########################################################################
# Some test problems
########################################################################

test: $(LIB)
	$(CC) $(CFLAGS) -c $(DRIVERS)/test.c -o $(DRIVERS)/test.o
	$(CC) $(LOPT) $(DRIVERS)/test.o $(CLFLAGS) -o test.ex

diff: $(LIB)
	$(CC) $(CFLAGS) -c $(DRIVERS)/ReactionAdvectionDiffusion.c -o $(DRIVERS)/ReactionAdvectionDiffusion.o
	$(CC) $(LOPT) $(DRIVERS)/ReactionAdvectionDiffusion.o $(CLFLAGS) -o diff.ex
########################################################################
# Clean up
########################################################################

.PHONY : clean distclean help

clean:
	rm -f $(CSRCDIR)/assemble/*.o
	rm -f $(CSRCDIR)/fem/*.o
	rm -f $(CSRCDIR)/grid/*.o
	rm -f $(CSRCDIR)/solver/*.o
	rm -f $(CSRCDIR)/utilities/*.o
	rm -f drivers/*.o

distclean:
	make clean
	rm -f lib/*.a
	rm -f include/*~
	rm -f *~ *.ex *.out
	rm -f $(CSRCDIR)/assemble/*~
	rm -f $(CSRCDIR)/fem/*~
	rm -f $(CSRCDIR)/grid/*~
	rm -f $(CSRCDIR)/solver/*~
	rm -f $(CSRCDIR)/utilities/*~

help:
	@echo "======================================================"
	@echo " 		   HAZMAT		             "
	@echo "======================================================"
	@echo " "
	@echo " make            : build all exe files "
	@echo " make headers    : build the header file automatically"
	@echo " make lib        : build library "
	@echo " make test	: build test.ex "
	@echo " make clean      : clean all obj files "
	@echo " make distclean  : clean all obj, exe, bak, out files "
	@echo " make help       : show this screen "
	@echo " "




#$(EXEeuler): $(OBJeuler) $(OBJS)
#	$(CC) $(LDFLAGS) $(OBJeuler) $(OBJS) $(LIBS) -o $(EXEeuler)

#$(EXEned): $(OBJned) $(OBJS)
#	$(FC) $(LDFLAGS) $(OBJned) $(OBJS) $(LIBS) -o $(EXEned) $(PROFFLAG)

#$(EXElap): $(OBJlap) $(OBJS)
#	$(FC) $(LDFLAGS) $(OBJlap) $(OBJS) $(LIBS) -o $(EXElap) $(PROFFLAG)

#$(EXEheat): $(OBJheat) $(OBJS)
#	$(FC) $(LDFLAGS) $(OBJheat) $(OBJS) $(LIBS) -o $(EXEheat) $(PROFFLAG)

#$(EXEstokes): $(OBJstokes) $(OBJS)
#	$(FC) $(LDFLAGS) $(OBJstokes) $(OBJS) $(LIBS) -o $(EXEstokes) $(PROFFLAG)

#$(EXEns): $(OBJns) $(OBJS)
#	$(FC) $(LDFLAGS) $(OBJns) $(OBJS) $(LIBS) -o $(EXEns) $(PROFFLAG)

#$(EXEmhd): $(OBJmhd) $(OBJS)
#	$(FC) $(LDFLAGS) $(OBJmhd) $(OBJS) $(LIBS) -o $(EXEmhd) $(PROFFLAG)

#$(EXEmoc): $(OBJmoc) $(OBJS)
#	$(FC) $(LDFLAGS) $(OBJmoc) $(OBJS) $(LIBS) -o $(EXEmoc)  $(PROFFLAG) 

#$(EXErt): $(OBJrt) $(OBJS)
#	$(FC) $(LDFLAGS) $(OBJrt) $(OBJS) $(LIBS) -o $(EXErt) $(PROFLAG)


#$(EXEfosls): $(OBJfosls) $(OBJS)
#	$(FC) $(LDFLAGS) $(OBJfosls) $(OBJS) $(LIBS) -o $(EXEfosls) $(PROFFLAG)

#$(EXEfoslsc): $(OBJfoslsc) $(OBJS)
#	$(FC) $(LDFLAGS) $(OBJfoslsc) $(OBJS) $(LIBS) -o $(EXEfoslsc) $(PROFFLAG)

#$(EXEmax): $(OBJmax) $(OBJS)
#	$(FC) $(LDFLAGS) $(OBJmax) $(OBJS) $(LIBS) -o $(EXEmax) $(PROFFLAG)

#$(EXEbdm): $(OBJbdm) $(OBJS)
#	$(FC) $(LDFLAGS) $(OBJbdm) $(OBJS) $(LIBS) -o $(EXEbdm) $(PROFFLAG)

#fheaders:	
#	/bin/cat $(SOURCES) \
#	| awk -f mkheaders.awk > $(MYINCLUDE)/fem_mhd_functs.h

#forts.h:	$(OBJF)
#		sed -f createh.sed $(SRC)/*.f > $(MYINCLUDE)/forts.h


#%.o:	%.c
#	$(CC) $(INCLUDE) $(CFLAGS) -o $@ -c  $< $(PROFFLAG)

#%.o:	%.f
#	$(FC) $(INCLUDE) $(FFLAGS) -o $@ -c  $< $(PROFFLAG)

#clean:
#	rm -f $(EXEtest) $(EXEeuler) $(EXEned) $(EXElap) $(EXEheat) $(EXEstokes) $(EXEns) $(EXEmhd) $(EXEmoc) $(EXErt) $(EXEfosls) $(EXEfoslsc) $(EXEmax) $(EXEbdm) $(EXEdg) $(SRC)/*.o $(DRIVERS)/*.o


#####################################################
# Last Modified 2012-07-21 --jha
####################################################

# Different Executable Programs
EXEtest = test
EXEeuler = euler
EXEned = NED
EXElap = Lap
EXEheat = Heat
EXEstokes = Stokes
EXEns = NS
EXEmhd = MHD
EXEmoc = MOC
EXErt = RT
EXEfosls = FOSLS
EXEfoslsc = CFOSLS
EXEmax = MAX
EXEbdm = BDM

# Generic Directories
DIR0 = .
MYINCLUDE = $(DIR0)/include
SRC = $(DIR0)/src
LIB = $(DIR0)/lib
DRIVERS = $(DIR0)/drivers

# Machine Specific Compilers and Libraries
# LINUX MACHINES
CC = gcc44
FC = gfortran44
CXX = g++44
ARPACK0 = $(LIB)/libarpack_Linux_64.a       #for Linux

# APPLE MACHINES
CC = gcc-mp-4.9
FC = gfortran-mp-4.9
CXX = g++-mp-4.9
CFLAGS = -g  
FFLAGS = -g -Wall -fno-second-underscore
LDFLAGS = -g -Wall
ARPACK0 = $(LIB)/libarpack_Darwin_x86_64.a #for MAC
# IF UMFPACK Needed include the following libraries and headers
SSDIR = /opt/local
MYDIR = /Users/jadler05/Documents/James/Research/CODE
# The METIS and HYPRE libraries (needed for the parallel examples)
METIS_DIR = $(MYDIR)/MFEM-SVN/metis-4.0
METISLIB = -L$(METIS_DIR) -lmetis
SSINCLUDE = -I$(SSDIR)/include
SSLIB = -L$(SSDIR)/lib -lsuitesparseconfig -lcholmod -lamd -lcolamd -lccolamd -lcamd -lspqr -lumfpack -lamd -lcxsparse
INCLUDE = -I$(MYINCLUDE) $(SSINCLUDE)
# Include all Libraries including UMFPACK, its dependencies, BLAS, and LAPACK
LIBS = $(SSLIB) $(ARPACK0) -lm -lblas -llapack -lgfortran

# c source and object files for drivers
SRCtest = $(DRIVERS)/7by7.c
OBJtest = $(DRIVERS)/7by7.o
SRCeuler = $(DRIVERS)/Euler.c
OBJeuler = $(DRIVERS)/Euler.o
SRCned = $(DRIVERS)/Nedmain.c
OBJned = $(DRIVERS)/Nedmain.o
SRClap = $(DRIVERS)/Laptest.c
OBJlap = $(DRIVERS)/Laptest.o
SRCheat = $(DRIVERS)/Heat.c
OBJheat = $(DRIVERS)/Heat.o
SRCstokes = $(DRIVERS)/Stokes.c
OBJstokes = $(DRIVERS)/Stokes.o
SRCns = $(DRIVERS)/NSmain.c
OBJns = $(DRIVERS)/NSmain.o
SRCmhd = $(DRIVERS)/MHDmain.c
OBJmhd = $(DRIVERS)/MHDmain.o
SRCmoc = $(DRIVERS)/NS_MOC.c
OBJmoc = $(DRIVERS)/NS_MOC.o
SRCrt = $(DRIVERS)/RTmain.c
OBJrt = $(DRIVERS)/RTmain.o
SRCfosls = $(DRIVERS)/Fosls.c
OBJfosls = $(DRIVERS)/Fosls.o
SRCfoslsc = $(DRIVERS)/Fosls_C.c
OBJfoslsc = $(DRIVERS)/Fosls_C.o
SRCmax = $(DRIVERS)/Maxwell_XJ.c
OBJmax = $(DRIVERS)/Maxwell_XJ.o
#SRCmax = $(DRIVERS)/Maxwell_short.c
#OBJmax = $(DRIVERS)/Maxwell_short.o
SRCbdm = $(DRIVERS)/BDMStokes.c
OBJbdm = $(DRIVERS)/BDMStokes.o

# Source and Object files for code
SOURCES = $(wildcard $(SRC)/*.c)
FSOURCES = $(wildcard $(SRC)/*.f)
OBJS=	$(SRC)/assemble.o \
	$(SRC)/assemble_short.o \
	$(SRC)/aux.o \
	$(SRC)/basis.o \
	$(SRC)/boundary.o \
	$(SRC)/data_input.o \
	$(SRC)/dumpdata.o \
	$(SRC)/fem_mappings.o \
	$(SRC)/gridinput.o \
	$(SRC)/quadrature.o \
	$(SRC)/solver.o \
	$(SRC)/timestep.o \
	$(SRC)/tri_stats.o \
	$(SRC)/callarpack.o \
	$(SRC)/eigsdrv1.o \
	$(SRC)/eigsdrv3.o \
	$(SRC)/spconvert0.o \
	$(SRC)/sgauss0.o

.PHONY: all

help:
	@echo -n $(SOURCES) $(FSOURCES) $(OBJS)

all: $(EXEmax)

$(EXEtest): $(OBJtest) $(OBJS)
	$(FC) $(LDFLAGS) $(OBJStest) $(OBJS) $(LIBS) -o $(EXEtest)

$(EXEeuler): $(OBJeuler) $(OBJS)
	$(CC) $(LDFLAGS) $(OBJeuler) $(OBJS) $(LIBS) -o $(EXEeuler)

$(EXEned): $(OBJned) $(OBJS)
	$(FC) $(LDFLAGS) $(OBJned) $(OBJS) $(LIBS) -o $(EXEned) $(PROFFLAG)

$(EXElap): $(OBJlap) $(OBJS)
	$(FC) $(LDFLAGS) $(OBJlap) $(OBJS) $(LIBS) -o $(EXElap) $(PROFFLAG)

$(EXEheat): $(OBJheat) $(OBJS)
	$(FC) $(LDFLAGS) $(OBJheat) $(OBJS) $(LIBS) -o $(EXEheat) $(PROFFLAG)

$(EXEstokes): $(OBJstokes) $(OBJS)
	$(FC) $(LDFLAGS) $(OBJstokes) $(OBJS) $(LIBS) -o $(EXEstokes) $(PROFFLAG)

$(EXEns): $(OBJns) $(OBJS)
	$(FC) $(LDFLAGS) $(OBJns) $(OBJS) $(LIBS) -o $(EXEns) $(PROFFLAG)

$(EXEmhd): $(OBJmhd) $(OBJS)
	$(FC) $(LDFLAGS) $(OBJmhd) $(OBJS) $(LIBS) -o $(EXEmhd) $(PROFFLAG)

$(EXEmoc): $(OBJmoc) $(OBJS)
	$(FC) $(LDFLAGS) $(OBJmoc) $(OBJS) $(LIBS) -o $(EXEmoc)  $(PROFFLAG) 

$(EXErt): $(OBJrt) $(OBJS)
	$(FC) $(LDFLAGS) $(OBJrt) $(OBJS) $(LIBS) -o $(EXErt) $(PROFLAG)


$(EXEfosls): $(OBJfosls) $(OBJS)
	$(FC) $(LDFLAGS) $(OBJfosls) $(OBJS) $(LIBS) -o $(EXEfosls) $(PROFFLAG)

$(EXEfoslsc): $(OBJfoslsc) $(OBJS)
	$(FC) $(LDFLAGS) $(OBJfoslsc) $(OBJS) $(LIBS) -o $(EXEfoslsc) $(PROFFLAG)

$(EXEmax): $(OBJmax) $(OBJS)
	$(FC) $(LDFLAGS) $(OBJmax) $(OBJS) $(LIBS) -o $(EXEmax) $(PROFFLAG)

$(EXEbdm): $(OBJbdm) $(OBJS)
	$(FC) $(LDFLAGS) $(OBJbdm) $(OBJS) $(LIBS) -o $(EXEbdm) $(PROFFLAG)

fheaders:	
	/bin/cat $(SOURCES) \
	| awk -f mkheaders.awk > $(MYINCLUDE)/fem_mhd_functs.h

forts.h:	$(OBJF)
		sed -f createh.sed $(SRC)/*.f > $(MYINCLUDE)/forts.h


%.o:	%.c
	$(CC) $(INCLUDE) $(CFLAGS) -o $@ -c  $< $(PROFFLAG)

%.o:	%.f
	$(FC) $(INCLUDE) $(FFLAGS) -o $@ -c  $< $(PROFFLAG)

clean:
	rm -f $(EXEtest) $(EXEeuler) $(EXEned) $(EXElap) $(EXEheat) $(EXEstokes) $(EXEns) $(EXEmhd) $(EXEmoc) $(EXErt) $(EXEfosls) $(EXEfoslsc) $(EXEmax) $(EXEbdm) $(EXEdg) $(SRC)/*.o $(DRIVERS)/*.o


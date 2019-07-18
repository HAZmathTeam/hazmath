#####################################################
# Last Modified 2017-12-17 --ltz
####################################################

HAZDIR = ../..
SRCDIR = .
##
OLDGRID = oldgridmake
GRID = gridmake

INCLUDE = .
LIBS = 
# Extension for the Executable Programs
EXTENSION = ex
# Machine Specific Compilers and Libraries
CC = gcc
FC = gfortran
CXX = g++
CFLAGS = -g -Wall
FFLAGS = -g -fno-second-underscore
ExtraFLAGS =

HAZLIB = $(HAZDIR)/lib

HAZINCLUDE = $(HAZDIR)/include

LIBS += -L$(HAZLIB) -lhazmath -lm -lgfortran #-lblas -llapack -lgfortran

HEADERS = $(INCLUDE)/grid_defs.h \
	$(INCLUDE)/grid_params.h

EXE = $(GRID).$(EXTENSION)
EXEold = $(OLDGRID).$(EXTENSION)
EXEin = ingrid.$(EXTENSION)
EXElex = lex.$(EXTENSION)

# Source and Object Files
OBJ =	$(SRCDIR)/main_grid.o	\
	$(SRCDIR)/newgrids.o \
	$(SRCDIR)/bitwise.o \
	$(SRCDIR)/new_refining.o		

OBJold =	$(SRCDIR)/$(GRID).o \
		$(SRCDIR)/new_refining.o \
		$(SRCDIR)/new_marking.o \
		$(SRCDIR)/features.o
##	$(SRCDIR)/auxfuncts.o 
##	$(SRCDIR)/io_geo.o \


OBJinp = main_input.o	sorting1.o input_grid1.o

OBJlex = lex1.o

.PHONY:  all

all: $(EXE) $(EXEold) $(EXEin) $(EXElex)

$(EXE):		$(OBJ)
		$(CC) $(LDFLAGS)  $(OBJ) -o $@ $(LIBS)

$(EXEold):	$(OBJold)
		$(CC) $(LDFLAGS)  $(OBJold) -o $@ $(LIBS)

$(EXEin):	$(OBJinp)
		$(CC) $(LDFLAGS)  $(OBJinp) -o $@

$(EXElex):	$(OBJlex)
		$(CC) $(LDFLAGS)  $(OBJlex) -o $@ $(LIBS)

%.o:	%.c	$(HEADERS)
	+$(CC) -I$(INCLUDE) -I$(HAZINCLUDE) $(CFLAGS) -o $@ -c $<

clean:
	+rm -rf *.ex $(SRCDIR)/*.o *.mod ./*.dSYM  
#########################################################

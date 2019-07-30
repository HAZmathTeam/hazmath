#####################################################
# Last Modified 2017-12-17 --ltz
####################################################

HAZDIR = ../..
SRCDIR = .
##
OLDGRID = oldgridmake
GRID = grid_generate
AMR = amr

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
EXEamr = $(AMR).$(EXTENSION)
EXEold = $(OLDGRID).$(EXTENSION)

# Source and Object Files
OBJ =	$(SRCDIR)/$(GRID).o
### 	$(SRCDIR)/macroelements.o	$(SRCDIR)/input_grid.o	$(SRCDIR)/amr_coords_z.o 

OBJamr =	$(SRCDIR)/$(AMR).o \
		$(SRCDIR)/refining.o \
		$(SRCDIR)/marking.o 
###		$(SRCDIR)/features.o

OBJold =	$(SRCDIR)/$(OLDGRID).o \
		$(SRCDIR)/new_refining.o \
		$(SRCDIR)/new_marking.o \
		$(SRCDIR)/features.o
.PHONY:  all

all: $(EXE) $(EXEold)	$(EXEamr)

$(EXE):		$(OBJ)
		$(CC) $(LDFLAGS)  $(OBJ) -o $@ $(LIBS)

$(EXEamr):	$(OBJamr)
		$(CC) $(LDFLAGS)  $(OBJamr) -o $@ $(LIBS)

$(EXEold):	$(OBJold)
		$(CC) $(LDFLAGS)  $(OBJold) -o $@ $(LIBS)

%.o:	%.c	$(HEADERS)
	+$(CC) -I$(INCLUDE) -I$(HAZINCLUDE) $(CFLAGS) -o $@ -c $<

clean:
	+rm -rf *.ex $(SRCDIR)/*.o *.mod ./*.dSYM  
#########################################################

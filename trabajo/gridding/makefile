#####################################################
# Last Modified 2017-12-17 --ltz
####################################################

HAZDIR = ../..
SRCDIR = .
##
OLDGRID = oldgridmake
GRID = grid_generate

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

# Source and Object Files
OBJ =	$(SRCDIR)/$(GRID).o \
	$(SRCDIR)/macroelements.o

OBJold =	$(SRCDIR)/$(OLDGRID).o \
		$(SRCDIR)/new_refining.o \
		$(SRCDIR)/new_marking.o \
		$(SRCDIR)/features.o
.PHONY:  all

all: $(EXE) $(EXEold)	

$(EXE):		$(OBJ)
		$(CC) $(LDFLAGS)  $(OBJ) -o $@ $(LIBS)

$(EXEold):	$(OBJold)
		$(CC) $(LDFLAGS)  $(OBJold) -o $@ $(LIBS)

%.o:	%.c	$(HEADERS)
	+$(CC) -I$(INCLUDE) -I$(HAZINCLUDE) $(CFLAGS) -o $@ -c $<

clean:
	+rm -rf *.ex $(SRCDIR)/*.o *.mod ./*.dSYM  
#########################################################

#####################################################
# Last Modified 2017-12-17 --ltz
####################################################

HAZDIR = ../..
SRCDIR = .
##
GRID = gridmake

INCLUDE = .
LIBS = 
# Extension for the Executable Programs
EXTENSION = ex
# Machine Specific Compilers and Libraries
CC = gcc
FC = gfortran
CXX = g++
CFLAGS = -g
FFLAGS = -g -fno-second-underscore
ExtraFLAGS =

HAZLIB = $(HAZDIR)/lib

HAZINCLUDE = $(HAZDIR)/include

LIBS += -L$(HAZLIB) -lhazmath -lm -lblas -llapack -lgfortran

HEADERS = $(INCLUDE)/grid_defs.h \
	$(INCLUDE)/grid_params.h

EXE = $(GRID).$(EXTENSION)

# Source and Object Files
OBJS = 	$(SRCDIR)/$(GRID).o \
	$(SRCDIR)/new_refining.o \
	$(SRCDIR)/new_marking.o \
	$(SRCDIR)/auxfuncts.o \
	$(SRCDIR)/io_geo.o \
	$(SRCDIR)/features.o

.PHONY:  all

all: $(EXE)

$(EXE):	$(OBJS)
	$(CC) $(LDFLAGS)  $(OBJS) -o $@ $(LIBS)


%.o:	%.c	$(HEADERS)
	+$(CC) -I$(INCLUDE) -I$(HAZINCLUDE) $(CFLAGS) -o $@ -c $<

clean:
	+rm -rf $(EXE) $(SRCDIR)/*.o *.mod ./*.dSYM  $(EXTRA_DEL)
#########################################################

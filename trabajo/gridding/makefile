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

LIBS += -L$(HAZLIB) -lhazmath -lm -lgfortran #-lblas -llapack -lgfortran

HEADERS = $(INCLUDE)/grid_defs.h \
	$(INCLUDE)/grid_params.h

EXE = $(GRID).$(EXTENSION)
EXEb = b.$(EXTENSION)
EXEin = ingrid.$(EXTENSION)
EXElex = lex.$(EXTENSION)

# Source and Object Files
OBJS = 	$(SRCDIR)/$(GRID).o \
	$(SRCDIR)/new_refining.o \
	$(SRCDIR)/new_marking.o \
	$(SRCDIR)/auxfuncts.o \
	$(SRCDIR)/io_geo.o \
	$(SRCDIR)/features.o

OBJb = 	newgrids.o \
	bitwise.o \
	new_refining.o

OBJinp = main_input.o

OBJlex = lex1.o

.PHONY:  all

all: $(EXE) $(EXEb) $(EXEin) $(EXElex)

$(EXE):		$(OBJS)
		$(CC) $(LDFLAGS)  $(OBJS) -o $@ $(LIBS)

$(EXEb):	$(OBJb)
		$(CC) $(LDFLAGS)  $(OBJb) -o $@ $(LIBS)

$(EXEin):	$(OBJinp)
		$(CC) $(LDFLAGS)  $(OBJinp) -o $@ $(LIBS)

$(EXElex):	$(OBJlex)
		$(CC) $(LDFLAGS)  $(OBJlex) -o $@ $(LIBS)

%.o:	%.c	$(HEADERS)
	+$(CC) -I$(INCLUDE) -I$(HAZINCLUDE) $(CFLAGS) -o $@ -c $<

clean:
	+rm -rf *.ex $(SRCDIR)/*.o *.mod ./*.dSYM  
#########################################################

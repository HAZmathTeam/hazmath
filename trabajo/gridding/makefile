#####################################################
# Last Modified 2017-12-17 --ltz
####################################################

HAZDIR = ../..
SRCDIR = .
##
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

LIBS += -L$(HAZLIB) -lhazmath -lm #-lblas -llapack -lgfortran

HEADERS = $(INCLUDE)/solve_estimate_mark.h

EXEamr = $(AMR).$(EXTENSION)

# Source and Object Files
OBJamr =	$(SRCDIR)/$(AMR).o
###		$(SRCDIR)/features.o

.PHONY:  all

all: $(EXE) $(EXEold)	$(EXEamr)

$(EXEamr):	$(OBJamr)
		$(CC) $(LDFLAGS)  $(OBJamr) -o $@ $(LIBS)

$(EXEold):	$(OBJold)
		$(CC) $(LDFLAGS)  $(OBJold) -o $@ $(LIBS)

%.o:	%.c	$(HEADERS)
	+$(CC) -I$(INCLUDE) -I$(HAZINCLUDE) $(CFLAGS) -o $@ -c $<

clean:
	+rm -rf *.ex $(SRCDIR)/*.o *.mod ./*.dSYM  
#########################################################

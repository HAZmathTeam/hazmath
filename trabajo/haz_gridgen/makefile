#####################################################
# Last Modified 2020-03-18 --ltz
####################################################
# SRC File name
EXE = ch.ex

# Machine Specific Compilers and Libraries
CC = gcc
CXX = g++
CFLAGS +=
FFLAGS += -fno-second-underscore

##################### no change should be needed below. ###########
# HAZMATH LIB and INCLUDE
HAZDIR = $(realpath ../../..)

HAZLIB = -L$(HAZDIR)/lib

INCLUDE = -I$(HAZDIR)/include

LIBS =  -lhazmath -lm -lblas -llapack

HEADERS =

HEADERS =

# Source and Object Files
OBJS =	convex_hull_n.o \
				visible_faces.o	\
				sort0.o

.PHONY:  all

all: $(EXE)

$(EXE):	$(OBJS)
	+$(CC) $(HAZLIB) $(OBJS) -o $@  $(LIBS)

%.o:	%.c
	+$(CC) $(INCLUDE) $(CFLAGS) -o $@ -c $<

clean:
	+rm -rf $(EXE) *.o *.mod output/* ./*.dSYM

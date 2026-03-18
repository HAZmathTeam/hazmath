#####################################################
# Liquid Crystals Pure Elastic example
####################################################

# SRC File name
SRCFILE = lcelastic
HEADERS =

INCLUDE =

LIBS =

# Switch between debug or optimized
CFLAGS = -g

## If an external libraries support is compiled in the hazmath library,
## uncomment  the corresponding lines below:
WITH_SUITESPARSE=1
WITH_LAPACK=1

sinclude ../../examples.mk

# Override HAZDIR since we are one level deeper than usual
override HAZDIR = $(realpath ../../..)
HAZLIB = -L$(HAZDIR)/lib -lhazmath
INCLUDE += -I$(HAZDIR)/include

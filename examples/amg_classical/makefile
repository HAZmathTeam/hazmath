#####################################################
# Classical Ruge-Stuben AMG + Incomplete Cholesky
# Last Modified 2026-03-04
####################################################

# Main program
SRCFILE = test_amg
EXTRA_SRCFILE = test_ilu_ichol

# Additional source files for the AMG library
AMG_SRCS = amg_wrapper.c

OBJS = $(AMG_SRCS:.c=.o)

HEADERS = amg_helping.h

INCLUDE =

LIBS =

# Switch between debug or optimized
CFLAGS = -O3
FFLAGS =

WITH_SUITESPARSE=1
WITH_LAPACK=1

sinclude ../examples.mk

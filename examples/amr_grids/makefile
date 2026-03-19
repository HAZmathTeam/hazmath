#####################################################
# Last Modified 2026-03-19
####################################################

# SRC File name
SRCFILE = amr_grids
HEADERS = supporting_$(SRCFILE).h

EXTRA_SRCFILE = test_dgs_adaptive

INCLUDE =

LIBS =

# Switch between debug or optimized
CFLAGS = -g


## If an external libraries support is compiled in the hazmath library,
## uncomment  the corresponding lines below:
WITH_SUITESPARSE=1
WITH_LAPACK=1

sinclude ../examples.mk

# Additional target: test_msh_roundtrip
test_msh_roundtrip.$(EXTENSION): test_msh_roundtrip.o
	+$(LINKER) $(CFLAGS) $(ExtraFLAGS) test_msh_roundtrip.o $(RPATH) -o $@ $(LIBS)

all: test_msh_roundtrip.$(EXTENSION)

clean: clean_roundtrip
clean_roundtrip:
	rm -f test_msh_roundtrip.$(EXTENSION) test_msh_roundtrip.o

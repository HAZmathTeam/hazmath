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

# Additional target: test_dgs_geometric
test_dgs_geometric.$(EXTENSION): test_dgs_geometric.o
	+$(LINKER) $(CFLAGS) $(ExtraFLAGS) test_dgs_geometric.o $(RPATH) -o $@ $(LIBS)

all: test_dgs_geometric.$(EXTENSION)

clean: clean_roundtrip
clean_roundtrip:
	rm -f test_msh_roundtrip.$(EXTENSION) test_msh_roundtrip.o
	rm -f test_dgs_geometric.$(EXTENSION) test_dgs_geometric.o

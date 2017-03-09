#####################################################
# Last Modified 2017-03-06 --ltz
####################################################

# SRC File name
SRCFILE = multigraph_solve

#Multigraph
MULTIGRAPH_DIR = $(realpath /Users/ltz1/packages/multigraph_2.0)
##########################################################
MULTIGRAPH_SRCDIR = $(MULTIGRAPH_DIR)/src
MULTIGRAPH_PRECISION_MOD = $(MULTIGRAPH_SRCDIR)/mthdef.mod
MULTIGRAPH_PRECISION_SRC = $(MULTIGRAPH_SRCDIR)/mg0.f

HEADERS_ADD = $(SRCFILE).h

OBJS_ADD = solver.o

EXTRA_DEL = $(MULTIGRAPH_SRCDIR)/mthdef.mod

sinclude ../examples.mk

%.o:	$(MULTIGRAPH_SRCDIR)/%.f $(MULTIGRAPH_PRECISION_MOD)
	+$(FC) $(INCLUDE) $(FFLAGS) -o $@ -c  $< 

$(MULTIGRAPH_PRECISION_MOD):	$(MULTIGRAPH_PRECISION_SRC)
		+$(FC)  $(CFLAGS) -c $<  $(INCLUDE) -o $@ 

# Generate an error message if the HAZmath library is not built
$(HAZLIB):
	$(error The HAZmath library is not built or is not up to date)


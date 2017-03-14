#####################################################
# Last Modified 2017-03-06 --ltz
####################################################

# SRC File name
SRCFILE = ConvectionDiffusion

HEADERS = $(SRCFILE).h

INCLUDE =

################comment the next two lines if you do not use mutigraph as solver.
WITH_MGRAPH = yes
MGRAPH_SRCDIR = $(realpath ../../../../multigraph_2.0/src)

# if you have compiled hazmath with suitesparse support, uncomment the
# two lines below.
#SSDIR = /usr/local
#LIBS = -L$(SSDIR)/lib -lsuitesparseconfig -lcholmod -lamd -lcolamd -lccolamd -lcamd -lspqr -lumfpack -lamd -lcxsparse 

#INCLUDE = -I$(SSDIR)/include

##########################################################

sinclude ../examples.mk

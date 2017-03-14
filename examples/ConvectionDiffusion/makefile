#####################################################
# Last Modified 2017-03-06 --ltz
####################################################

# SRC File name
SRCFILE = ConvectionDiffusion
#SRCFILE = multigraph_solve

# if you have compiled hazmath with suitesparse support, uncomment the
# two lines below.
SSDIR = /usr/local
SSLIB = -L$(SSDIR)/lib -lsuitesparseconfig -lcholmod -lamd -lcolamd -lccolamd -lcamd -lspqr -lumfpack -lamd -lcxsparse 

LIBS_ADD = $(SSLIB) 

##########################################################
HEADERS_ADD = $(SRCFILE).h

WITH_MGRAPH = yes

sinclude ../examples.mk

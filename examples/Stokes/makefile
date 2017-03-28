#####################################################
# Last Modified 2017-02-11 --jha
####################################################

# SRC File name
SRCFILE = Stokes

HEADERS = StokesData.h StokesSystem.h

INCLUDE =

# if you have compiled hazmath with suitesparse support, uncomment the
# two lines below.
SSDIR = /usr/local
LIBS = -L$(SSDIR)/lib -lsuitesparseconfig -lcholmod -lamd -lcolamd -lccolamd -lcamd -lspqr -lumfpack -lamd -lcxsparse 

INCLUDE = -I$(SSDIR)/include

sinclude ../examples.mk

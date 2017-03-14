#####################################################
# Last Modified 2017-02-11 --jha
####################################################

# SRC File name
SRCFILE = Stokes
HEADERS_DO = StokesData.h StokesSystem.h

# if you have compiled hazmath with suitesparse support, uncomment the
# two lines below.
SSDIR = /usr/local
SSLIB = -L$(SSDIR)/lib -lsuitesparseconfig -lcholmod -lamd -lcolamd -lccolamd -lcamd -lspqr -lumfpack -lamd -lcxsparse

LIBS_ADD = $(SSLIB)

sinclude ../examples.mk

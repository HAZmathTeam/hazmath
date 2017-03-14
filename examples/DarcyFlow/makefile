#####################################################
# Last Modified 2017-03-09 --ltz
####################################################

# SRC File name
SRCFILE = DarcyFlow
HEADERS_ADD = DarcySystem.h DarcyData.h

# if you have compiled hazmath with suitesparse support, uncomment the
# two lines below.
SSDIR = /usr/local
SSLIB = -L$(SSDIR)/lib -lsuitesparseconfig -lcholmod -lamd -lcolamd -lccolamd -lcamd -lspqr -lumfpack -lamd -lcxsparse

LIBS_ADD = $(SSLIB)

sinclude ../examples.mk

# Generate an error message if the HAZmath library is not built
$(HAZLIB):
	$(error The HAZmath library is not built or is not up to date)





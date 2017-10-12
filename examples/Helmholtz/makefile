#####################################################
# Last Modified 2017-03-08 --ltz
####################################################
# SRC File name
SRCFILE = helmholtz

HEADERS = 

INCLUDE =

# Switch between debug or optimized
CFLAGS = -g
FFLAGS = -g

# if you have compiled hazmath with suitesparse support, uncomment the
# two lines below.
WITH_SUITESPARSE=0

ifeq ($(WITH_SUITESPARSE),1)
	CFLAGS += -DWITH_SUITESPARSE=1
	SSDIR = /opt/local
	LIBS = -L$(SSDIR)/lib -lsuitesparseconfig -lcholmod -lamd -lcolamd -lccolamd -lcamd -lspqr -lumfpack -lamd -lcxsparse 

	INCLUDE = -I$(SSDIR)/include
endif

sinclude ../examples.mk



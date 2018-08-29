#####################################################
# Last Modified 2017-04-03 --jha & xh
####################################################
# SRC File name
SRCFILE = haz_solve

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
	SSDIR = /usr/local
	LIBS = -L$(SSDIR)/lib -lsuitesparseconfig -lcholmod -lamd -lcolamd -lccolamd -lcamd -lspqr -lumfpack -lamd -lcxsparse 

	INCLUDE = -I$(SSDIR)/include
endif

sinclude ../examples.mk

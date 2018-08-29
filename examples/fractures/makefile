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
WITH_SUITESPARSE=1
LIBS += -L$(HAZLIB) -lhazmath -lm -lblas -llapack 

# if you have compiled hazmath with suitesparse support, uncomment the
# two lines below.
WITH_SUITESPARSE=1
INCLUDESSP=
ifeq ($(WITH_SUITESPARSE),1)
	CFLAGS += -DWITH_SUITESPARSE=1
	SSDIR = /usr/lib/x86_64-linux-gnu
	INCLUDESSP = /usr/include/suitesparse
	LIBS += -lsuitesparseconfig -lcholmod -lamd -lcolamd -lccolamd -lcamd -lspqr -lumfpack -lamd -lcxsparse 
endif

sinclude ../examples.mk

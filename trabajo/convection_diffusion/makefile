#####################################################
# Last Modified 2017-03-06 --ltz
####################################################

# SRC File name
SRCFILE = convectiondiffusion

HEADERS = $(SRCFILE).h

OBJS += eafe.o

INCLUDE =

################comment the next two lines if you do not use mutigraph as solver.
WITH_MGRAPH = yes
MGRAPH_SRCDIR = $(realpath /home/ltz1/packages/multigraph_2.1)

# if you have compiled hazmath with suitesparse support, uncomment the
# two lines below.
#SSDIR = /usr/local
#LIBS = -L$(SSDIR)/lib -lsuitesparseconfig -lcholmod -lamd -lcolamd -lccolamd -lcamd -lspqr -lumfpack -lamd -lcxsparse 

#INCLUDE = -I$(SSDIR)/include

##########################################################
LIBS = -lgfortran

sinclude conv_diff.mk


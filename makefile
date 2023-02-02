#######################################################################
#
########################################################################
# TOP LEVEL HAZMATH Makefile: Calls cmake to configure and build
# the HAZMATH library and examples.
# 
#   Probably you will *NOT NEED TO CHANGE* this  top level Makefile.
#
#   USER DEFINED OPTIONS GO IN "hazmath.mk" which is included by this
#   makefile. Copy "hazmath.mk.example file to "hazmath.mk", edit it to
#   adjust the settings to your liking, and then type "make help" to
#   see how to configure/build HAZMATH.
# 
#  Modified   		2015-08-08   	--ltz
#  Added Doxygen	2017-03-12	--Xiaozhe Hu
#  Modified cflags to suppress warnings under linux 2017-03-12 --ltz
########################################################################
sinclude haz_config/hazmath.mk

cflags="-Wno-unused-result"
cxxlags="-Wno-unused-result"
#fflags="-Wno-unused-result" this is not valid in fortran
ifeq ($(debug),yes)
	cflags="-Wall -g"
	cxxflags="-Wall -g"
	fflags="-Wall -g"
endif
#
####################  User Changes UP TO HERE   ########################

# Let cmake do the configuration. Set up a build dir
build_dir=BUILD_HAZ

CONFIG_FLAGS=-DCMAKE_RULE_MESSAGES=ON

ifeq ($(verbose),yes)
    CONFIG_FLAGS+=-DCMAKE_VERBOSE_MAKEFILE=ON
else
    CONFIG_FLAGS+=-DCMAKE_VERBOSE_MAKEFILE=OFF
endif

ifeq ($(debug),yes)
    CONFIG_FLAGS+=-DCMAKE_BUILD_TYPE=DEBUG
else
    CONFIG_FLAGS+=-DCMAKE_BUILD_TYPE=RELEASE
endif

ifeq ($(shared),yes)
    CONFIG_FLAGS+=-DSHARED=$(shared)
endif

ifeq ($(openmp),yes)
    CONFIG_FLAGS+=-DUSE_OPENMP=$(openmp)
endif

ifeq ($(blas), yes)
    CONFIG_FLAGS+=-DUSE_BLAS=$(blas)
    CONFIG_FLAGS+=-DBLAS_DIR=$(blas_dir)
endif

ifeq ($(lapack), yes)
    CONFIG_FLAGS+=-DUSE_LAPACK=$(lapack)
    CONFIG_FLAGS+=-DLAPACK_DIR=$(lapack_dir)
endif


ifeq ($(suitesparse), yes)
    CONFIG_FLAGS+=-DUSE_SUITESPARSE=$(suitesparse)
    CONFIG_FLAGS+=-DSUITESPARSE_DIR=$(suitesparse_dir)
endif

ifeq ($(hdf5), yes)
    CONFIG_FLAGS+=-DUSE_HDF5=$(hdf5)
    CONFIG_FLAGS+=-DHDF5_DIR=$(hdf5_dir)
endif

ifeq ($(fftw3), yes)
    CONFIG_FLAGS+=-DUSE_FFTW3=$(fftw3)
    CONFIG_FLAGS+=-DFFTW3_DIR=$(fftw3_dir)
endif


ifeq ($(haznics), yes)
    CONFIG_FLAGS+=-DUSE_HAZNICS=$(haznics)
#    CONFIG_FLAGS+=-DHAZNICS_DIR=$(haznics_dir)
endif

ifeq ($(doxygen),yes)
    CONFIG_FLAGS+=-DUSE_DOXYGEN=$(doxygen)
endif

ifeq ($(multigraph), yes)
    CONFIG_FLAGS+=-DUSE_MULTIGRAPH=$(multigraph)
    CONFIG_FLAGS+=-DMULTIGRAPH_DIR=$(multigraph_dir)
endif

ifeq ($(matlab), yes)
    CONFIG_FLAGS+=-DUSE_MATLAB=$(matlab)
endif

CONFIG_FLAGS+=-DADD_CFLAGS=$(cflags)
CONFIG_FLAGS+=-DADD_CXXFLAGS=$(cxxflags)
CONFIG_FLAGS+=-DADD_FFLAGS=$(fflags)

all clean headers docs:
	@if [ ! -f $(build_dir)/Makefile ] ; then \
		echo "Configuration not found! Please perform configuration first."; \
		echo "See the following help screen for usage ..."; \
		echo " "; \
		cat haz_docs/HELP.txt; \
	else \
	  	make -C $(build_dir) $@ ; \
	fi

install:	headers	
	@if [ ! -f $(build_dir)/Makefile ] ; then \
		echo "Configuration not found! Please perform configuration first."; \
		echo "See the following help screen for usage ..."; \
		echo " "; \
		cat haz_docs/HELP.txt; \
	else \
	  	make -C $(build_dir) install ; \
	fi
config: distclean
	@if [ ! -f ./haz_config/hazmath.mk ] ; then \
		echo "***ERROR: haz_config/hazmath.mk is missing...." ; \
		echo "   1. Copy \"haz_config/hazmath.mk.example\" to \"haz_config/hazmath.mk\"." ; \
		echo "   2. Adjust \"haz_config/hazmath.mk\" for your system or leave it with the default settings" ; \
		echo "   3. Run \"make config\" again." ; \
	else \
		mkdir -p $(build_dir) ; \
		cd $(build_dir) && cmake $(CURDIR) $(CONFIG_FLAGS) ; \
	fi
	@-echo " "
	@-echo "--------------------------------------------------------"
	@-echo "If SUCCESS, run 'make install' to install the library."
	@-echo "--------------------------------------------------------"
	@-echo " "

uninstall:
	@if [ ! -f $(build_dir)/install_manifest.txt ]; then \
		echo "Installation manifest not found! Nothing to uninstall."; \
		echo "See the following help screen for usage ..."; \
		echo " "; \
		cat haz_docs/HELP.txt; \
	else \
		xargs rm < $(build_dir)/install_manifest.txt; \
		rm -rf $(build_dir)/install_manifest.txt \
		       doc/htdocs; \
	fi

distclean:
	@-rm -rf $(build_dir)
	@-rm -rf ./lib/*	
	@-rm -rf ./swig_files/*	
	@-find . -name '*~' -exec rm {} \;

help:
	@clear
	@cat haz_docs/HELP.txt

.PHONY:	all config clean distclean install uninstall docs headers help


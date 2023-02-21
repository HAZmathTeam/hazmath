#######################################################################
#                Simple Finite Element Package (HAZMATH)
#
################# User Defined Configuration Options #################
#
# 1. Copy this file to a file named "hazmath.mk".
# 2. Edit "hazmath.mk" to adjust options/settings for your system
#    following the directions below.
# 3. Type "make help" to see all build and configuration options.
########################################################################
#
# The default setting for build type for HAZMATH is RELEASE. The RELEASE
# build type by default has the "-O3". You may adjust the optimization
# compilation options according to your hardware and software setting.
# For example, on a macbook pro with Intel i7, best options could be
# "-Ofast -march=corei7 -mtune=corei7".
#
# If you want to work with build type DEBUG, then uncomment the next
# line (to include "-Wall -g")
#
# debug=yes
#
# The default setting for vebosity level for HAZMATH is verbose=no. If you
# want to increase verbosity level, uncomment the next line:
#
# verbose=yes
#
# By default, HAZMATH generates static libraries. If you need to generate
# shared libs instead of static libs, uncomment the next line:
#
# shared=yes
#
# You may use multithread version after you enable OpenMP support. To
# setup the environment, you need
#  >> export OMP_NUM_THREADS=4 (for bash)
#  >> setenv OMP_NUM_THREADS 4 (for tcsh)
# If you want to compile with OpenMP support, uncomment the next line:
#
# openmp=yes
#
# These user options can also be applied as make command line options.
# For example, to enforce the debug compiling options:
#
# make config debug=yes
#
#-------------------------------------------------------------------------
#
# HAZMATH uses the command-line Doxygen to generate a reference manual.
# If you want to use the GUI of Doxgen instead of command-line,
# uncomment the next line:
#
# doxygen=yes
#
#-------------------------------------------------------------------------
# If you want to use the BLAS
#
# blas=yes
# blas_dir=/home/whatever/blas-xxx-yyy
#
# IMPORTANT: The BLAS libraries/(include files) are expected to be
# fouund either under the path specified above of in the system
# standard paths for libraries and header files.
# -------------------------------------------------------------------------
# If you want to use the LAPACK
#
# lapack=yes
# lapack_dir=/home/whatever/lapack-xxx-yyy
# FOR MAC
# lapack_dir=/opt/homebrew/opt/lapack
#
# IMPORTANT: The LAPACK libraries/(include files) are expected to be
# fouund either under the path specified above of in the system
# standard paths for libraries and header files.
# -------------------------------------------------------------------------
# If you want to use the SuiteSparse package, uncomment the next line
# (and read carefully the instructons below it):
#
# suitesparse=yes
#
# If you have installed SuiteSparse from source or for some other
# reason you want to specify the path to SuiteSparse libraries and
# header files, uncomment and edit the definition of "suitesparse_dir"
# below (and continue reading...)
#
# suitesparse_dir=/opt/local/lib
# suitesparse_dir=/home/xiaozhehu/Work/Projects/HAZMATH/SuiteSparse-5.3.0/SuiteSparse/lib
#
# IMPORTANT:
# This defines the path to the SuiteSparse library and include files.
# These are expected to be found in $(suitesparse_dir)/lib and
# $(suitesparse_dir)/include or in the system standard paths for libraries
# and header files.
# -------------------------------------------------------------------------
# If you want to use the HDF5 data library
#
# hdf5=yes
#
# -------------------------------------------------------------------------
# If you want to use the FFTW3 library
#
# fftw3=yes
#
# -------------------------------------------------------------------------
# If you want to use the interface with MATLAB, uncomment the next line:
#
# matlab=yes
#
# -------------------------------------------------------------------------
####################  User Defined Compiler Flags  #####################
#cflags="-funroll-loops -funswitch-loops"
#cxxflags="-funroll-loops -funswitch-loops"
#fflags="-funroll-loops -funswitch-loops"

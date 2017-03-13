#!/bin/sh
# 
# usage is tarmk.sh 
# Edit to include various directories and files into tarball that will be exhanged to users

VERSION="v00"
FILENAME="HAZMATH_$(date +%Y_%m_%d)_$VERSION.tar"

echo "Creating Files for HAZMATH $VERSION on $(date +%Y_%m_%d)"

cd ..

# CMAKE and CONFIG files
tar -cvf $FILENAME ./CMakeLists.txt \
                   ./haz_config/hazmath.mk \
                   ./makefile

# License and README
tar -rv -f $FILENAME ./LICENSE.txt \
                     ./README

# Docs
tar -rv -f $FILENAME ./haz_docs/HELP.txt \
                     ./haz_docs/hazmath.Doxygen.cnf.in

# Shell Utilities
tar -rv -f $FILENAME ./haz_shutils/headmk.sh \
                     ./haz_shutils/mkheaders.awk

# Includes
tar -rv -f $FILENAME ./include/*.h

# Library (just create folder)
tar --exclude=*.a --exclude=.git* -rv -f $FILENAME ./lib/

# Source Files
tar -rv -f $FILENAME ./src/assemble/*.c \
                     ./src/fem/*.c \
                     ./src/grid/*.c \
                     ./src/nonlinear/*.c \
                     ./src/solver/*.c \
                     ./src/solver/*.inl \
                     ./src/timestepping/*.c \
                     ./src/utilities/*.c \

# Examples
tar -rv -f $FILENAME ./examples/README \
                     ./examples/examples.mk \
                     ./examples/HDEquation/*.c ./examples/HDEquation/input.dat ./examples/HDEquation/makefile \
                     ./examples/HeatEquation/*.c ./examples/HeatEquation/input.dat ./examples/HeatEquation/makefile \
                     ./examples/Stokes/*.c ./examples/Stokes/*.h ./examples/Stokes/input.dat ./examples/Stokes/makefile \

# Grids
tar -rv -f $FILENAME ./examples/grids/1D/*.haz \
                     ./examples/grids/2D/unitSQ_hp*.haz \
                     ./examples/grids/3D/unitCUBE_v*.haz

mv $FILENAME ./haz_shutils/
cd ./haz_shutils
gzip -v $FILENAME

echo "$FILENAME.gz created."

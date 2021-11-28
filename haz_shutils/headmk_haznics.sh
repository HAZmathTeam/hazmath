#!/bin/sh
# 
# usage is headmk_swig.sh ${CMAKE_CURRENT_SOURCE_DIR} where we assume that
# ${CMAKE_CURRENT_SOURCE_DIR} has subdirectories src and include and
# the function names from src/haznics/*.c are put into the include/hazmath_swig.h
# commented as is not needed for now (ltz, 20190730) FORTRAN: the
#set +x

cat $1/src/haznics/*.c \
	| awk -v name="haznics_add.h" -f mkheaders_simple.awk > $1/include/haznics_add.h

#set -x

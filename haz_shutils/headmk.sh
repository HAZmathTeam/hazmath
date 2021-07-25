#!/bin/sh
# 
# usage is headmk.sh ${CMAKE_CURRENT_SOURCE_DIR} where we assume that
# ${CMAKE_CURRENT_SOURCE_DIR} has subdirectories src and include and
# the function names from src/*/*.c are put into the include/hazmath.h
# commented as is not needed for now (ltz, 20190730)
#set +x

cat $1/src/assemble/*.c $1/src/fem/*.c \
    $1/src/mesh/*.c $1/src/nonlinear/*.c \
    $1/src/solver/*.c  $1/src/utilities/*.c \
    $1/src/timestepping/*.c $1/src/interfaces/*.c  \
    $1/src/amr/*.c $1/src/graphs/*.c \
    $1/src/eigen/*.c \
    $1/src/approximation/*.c \
	| awk -v name="hazmath.h" -f mkheaders.awk > $1/include/hazmath.h

#cat $1/src/haznics/*.CC \
#	| awk -v name="helpder.hidden" -f mkheaders.awk >> $1/include/hazmath.h

#set -x

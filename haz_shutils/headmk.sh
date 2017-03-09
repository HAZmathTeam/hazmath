#!/bin/sh
# 
# usage is headmk.sh ${CMAKE_CURRENT_SOURCE_DIR} where we assume that
# ${CMAKE_CURRENT_SOURCE_DIR} has subdirectories src and include and
# the function names from src/*/*.c are put into the include/hazmath.h
set +x
/bin/cat $1/src/assemble/*.c $1/src/fem/*.c \
    $1/src/grid/*.c $1/src/nonlinear/*.c $1/src/solver/*.c  $1/src/utilities/*.c $1/src/timestepping/*.c $1/src/interfaces/*.c \
	| awk -v name="hazmath.h" -f mkheaders.awk > $1/include/hazmath.h
set -x

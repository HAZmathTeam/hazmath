#!/bin/bash

make clean; make
set -x
gcc -c -I../../include main_xd_1d.c ; gcc -L. -L../../lib -Wl,-rpath=../../lib -Wl,-rpath=./ -o main_xd_1d.ex main_xd_1d.o -lxd_1d_ -lhazmath
set +x

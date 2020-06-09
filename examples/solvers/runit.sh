#!/bin/sh
make clean ; make solvers.ex
for i in `ls -1 matrices/A_*` ; do
  j=`echo $i | sed -e 's/^matrices\/A_//' -e 's/\.txt$//'`;
  ofile="output/o_$j"
  afile="matrices/A_$j.txt"
  bfile="matrices/b_$j.txt"
  if [ -f $afile -a -f $bfile -a -d output ]; then
    set -x
    ./solvers.ex $afile $bfile > $ofile
    set +x
  fi
done

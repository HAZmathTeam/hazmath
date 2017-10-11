#!/bin/sh
theD="1D"

HAZHOME=..

gcc -O -Wall -Wno-unused-variable -o t_grid -L$HAZHOME/lib -I$HAZHOME/include trans_grids.c -lhazmath -lm


for i in `find $HAZHOME/examples/old-grids/$theD \( -name '*.grd' -a -type f \)`
do
    j=`echo $i|sed -e s/old-grids/grids/`
    ./t_grid $i $j
    diff $i $j
done

#for i in `find ../examples/old-grids/$theD \( -name *.haz -a -type f \)`
#do
#    j=`echo $i|sed -e s/old-grids/grids/`
#    echo $i $j
#done

#!/bin/sh
theD="1D"

HAZHOME=..

gcc -O -Wall -Wno-unused-variable -o t_grid -L$HAZHOME/lib -I$HAZHOME/include trans_grids.c -lhazmath -lm


#for i in `find $HAZHOME/examples/old-grids/$theD \( -name '*.grd' -a -type f \)`
#do
#    j=`echo $i|sed -e s/old-grids/grids/`
#    ./t_grid $i $j
    ./t_grid  $HAZHOME/examples/old-grids/3D/unitCUBE_n129.haz $HAZHOME/examples/grids/3D/unitCUBE_n129.haz
    diff $HAZHOME/examples/old-grids/3D/unitCUBE_n129.haz $HAZHOME/examples/grids/3D/unitCUBE_n129.haz
    ./t_grid  $HAZHOME/examples/old-grids/3D/unitCUBE_n65.haz $HAZHOME/examples/grids/3D/unitCUBE_n65.haz
    diff $HAZHOME/examples/old-grids/3D/unitCUBE_n65.haz $HAZHOME/examples/grids/3D/unitCUBE_n65.haz
#done

#for i in `find ../examples/old-grids/$theD \( -name *.haz -a -type f \)`
#do
#    j=`echo $i|sed -e s/old-grids/grids/`
#    echo $i $j
#done

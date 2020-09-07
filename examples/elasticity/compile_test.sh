set -x
gcc -o test_faces faces.c -I../../include -L../../lib -lm -lhazmath -lsuitesparseconfig -lcholmod -lamd -lcolamd -lccolamd -lcamd -lspqr -lumfpack -lamd -lcxsparse
set +x



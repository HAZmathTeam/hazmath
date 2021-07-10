gcc -g  -fPIC -c helpers.c -o helpers.o  
gcc -g   helpers.c   main.c -o app.exe

./app.exe


# swig stuff 
swig -python -py3 -I. helpers.i
gcc -g -fPIC -c helpers.c -o helpers.o 
gcc -g -fPIC -c pyhelpers.c -I/usr/include/python3.8 -I/usr/lib/python3/dist-packages/numpy/core/include/numpy  -o pyhelpers.o
gcc -fPIC -c helpers_wrap.c -I/usr/include/python3.8 -I/usr/lib/python3/dist-packages/numpy/core/include/numpy  -o helpers_wrap.o
ld -fPIC -shared helpers.o -L../../lib pyhelpers.o helpers_wrap.o -o _helpers.so -lpython3.8

python3 test.py 




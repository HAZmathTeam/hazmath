#!/bin/bash

# Update cmake if needed
echo "Installing cmake 3.20.2..."
wget https://github.com/Kitware/CMake/releases/download/v3.20.2/cmake-3.20.2-linux-x86_64.sh
chmod +x cmake-3.20.2-linux-x86_64.sh
yes | bash cmake-3.20.2-linux-x86_64.sh
# ln -s /link/to/cmake-3.20.2-linux-x86_64/bin/* /usr/local/bin
echo "... Done."

# Install swig4.0
echo "Installing swig 4.0.2..."
wget http://prdownloads.sourceforge.net/swig/swig-4.0.2.tar.gz
tar xvzf swig-4.0.2.tar.gz
cd swig-4.0.2
wget https://ftp.pcre.org/pub/pcre/pcre-8.44.tar.gz
./configure
make
make install
cd ..
echo "... Done."

# Install suitesparse if needed
# echo "Installing suitesparse"
# sudo apt-get update
# sudo apt install -y libsuitesparse-dev
# echo "Done installing suitesparse"

# Install fenics
echo "Installing fenics..."
sudo apt-get install software-properties-common
sudo add-apt-repository ppa:fenics-packages/fenics
sudo apt-get update
sudo apt-get install fenics
echo "... Done."

# Install pyamg -- i don't think this is needed
# git clone https://github.com/pyamg/pyamg.git
# cd pyamg
# python3 setup.py install --user
# cd ..

# Install quadpy (xii dependency - which avoids scipy.linalg import eigh_tridiagonal)
echo "Installing quadpy..."
git clone https://github.com/nschloe/quadpy.git
cd quadpy
git checkout v0.12.10
python3 setup.py install --user
cd ..
echo "... Done."

# Install fenics_ii (xii)
echo "Installing fenics_ii..."
git clone https://github.com/MiroK/fenics_ii.git
cd fenics_ii
git fetch --all
git checkout hsfrac-minimal
cd ..

export PYTHONPATH="./fenics_ii/":"$PYTHONPATH"
echo "... Done."

# Install hsmg
echo "Installing hsmg..."
git clone https://github.com/MiroK/hsmg.git
cd hsmg
git fetch --all
# python3 setup.py install --user
source setup.rc
cd ..
echo "... Done."

# export PYTHONPATH="./hsmg/":"$PYTHONPATH"

# Install cbc.block
echo "Installing cbc.block..."
git clone https://mirok-w-simula@bitbucket.org/mirok-w-simula/cbc.block.git
cd cbc.block
python3 setup.py install --user
cd ..
echo "... Done."

# Install ulfy
echo "Installing ulfy..."
git clone https://github.com/MiroK/ulfy.git
cd ulfy
git checkout python3
python3 setup.py install --user
cd ..
echo "... Done."

# Install hazmath (with the python interface haznics=yes)
echo "Installing hazmath (with haznics)..."
#git clone https://github.com/HAZmathTeam/hazmath.git
#cd hazmath
make -C ../.. config shared=yes suitesparse=yes lapack=yes haznics=yes
make -C ../.. install
##cd ..

export PYTHONPATH="../../swig_files/":"$PYTHONPATH"
echo "... Done."
# I don't know if below is necessary
# pip3 install --user networkx
# pip3 install --user gmsh

# sudo apt install --reinstall python*-decorator
# sudo apt install --reinstall python*-tabulate

#!/bin/bash
# A simple shell script to get haznics python dependencies
# Setup by Ana Budisa - 2021-07-15

# Install fenics
echo "Installing fenics..."
sudo apt-get install software-properties-common
sudo add-apt-repository ppa:fenics-packages/fenics
sudo apt-get update
sudo apt-get install fenics
echo "... Done."

# Install quadpy (xii dependency which avoids eigh_tridiagonal)
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
python3 setup.py install --user
cd ..
echo "... Done."

# Install hsmg
# echo "Installing hsmg..."
# git clone https://github.com/MiroK/hsmg.git
# cd hsmg
# git fetch --all
# source setup.rc
# cd ..
# echo "... Done."
## export PYTHONPATH="./hsmg/":"$PYTHONPATH"

# Install cbc.block
echo "Installing cbc.block..."
git clone https://bitbucket.org/fenics-apps/cbc.block.git
cd cbc.block
python3 setup.py install --user
cd ..
echo "... Done."

# Install ulfy
echo "Installing ulfy..."
git clone https://github.com/MiroK/ulfy.git
cd ulfy
python3 setup.py install --user
cd ..
echo "... Done."

# Install hazmath (with the python interface haznics=yes)
echo "Installing hazmath (with haznics)..."
make -C ../.. config shared=yes suitesparse=yes lapack=yes haznics=yes swig=yes
make -C ../.. install

echo "... Done."


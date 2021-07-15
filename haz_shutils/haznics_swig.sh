#!/bin/bash
# A simple shell script to get haznics dependencies
# Ana Budisa - 2021-07-15

# Update cmake if needed
echo "Installing cmake v >= 3.20..."
sudo snap install cmake
echo "... Done."

# Install swig4.0
echo "Installing swig 4.0.2..."
wget http://prdownloads.sourceforge.net/swig/swig-4.0.2.tar.gz
tar xvzf swig-4.0.2.tar.gz
cd swig-4.0.2
wget https://ftp.pcre.org/pub/pcre/pcre-8.44.tar.gz
./configure
make
sudo make install
cd ..
echo "... Done."

# Install suitesparse if needed
echo "Installing suitesparse..."
sudo apt-get update
sudo apt install -y libsuitesparse-dev
echo "... Done."
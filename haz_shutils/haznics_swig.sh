#!/bin/bash
# A simple shell script to get haznics dependencies
# Ana Budisa - 2021-07-15

# Update cmake if needed
echo "Installing cmake v >= 3.16..."
sudo apt-get update
sudo apt-get install cmake
echo "... Done."

# Install swig4.0
echo "Installing swig v >= 4.0..."
sudo apt-get install swig
echo "... Done."

# Install suitesparse if needed
echo "Installing suitesparse..."
sudo apt install -y libsuitesparse-dev
echo "... Done."

# Install python-dev
echo "Installing python3-dev..."
sudo apt-get install -y python3-dev
echo "... Done."

# Install numpy
echo "Installing numpy on python3..."
sudo python3 -m pip install --upgrade pip
sudo pip install numpy
echo "... Done."
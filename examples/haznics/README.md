### Requirements
- swig >= 4.0.0 (check [this](http://www.swig.org/download.html) for software download and [this](https://github.com/swig/swig/wiki/Getting-Started) for installing)
- python >= 3.5 (with numpy >= 1.13.3)
- cmake >= 3.12

### Install additional packages 
...required by the haznics demos: 
- Edit `setup.rc` and, if needed, `install_haznics.sh` to match your system
- Source the setup file (this executes install_haznics.sh and installs all package dependencies)
```
source setup.rc
```

### Running tests
```
cd tests/
python3 test_haznics.py
```

### Running demos, e.g. Poisson equation
```
python3 poisson.py
```

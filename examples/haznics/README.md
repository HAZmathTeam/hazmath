### Requirements
- swig >= 4.0.0 (check [this](http://www.swig.org/download.html) for software download and [this](https://github.com/swig/swig/wiki/Getting-Started) for installing)
- python >= 3.5 (with numpy >= 1.13.3)
- cmake >= 3.12

### Install additional packages (required by HAZniCS demos): 
- Edit, if needed, `install_haznics.sh` to match your system and execute with 
```
/bin/sh ./install_haznics.sh
``` 
(this installs all package dependencies: [FEniCS v2019.1.0](https://fenicsproject.org/download/archive/), [FEniCS_ii](https://github.com/MiroK/fenics_ii), [cbc.block](https://bitbucket.org/fenics-apps/cbc.block/) )
- add HAZniCS library to PYTHONPATH with 
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

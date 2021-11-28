Install additional packages required by the haznics demos: 
- edit setup.rc and, if needed, install_haznics.sh to match your system

### Source the setup file (this executes install_haznics.sh and installs all package dependencies)
```
source setup.rc
```

### Running tests:
```
cd tests/
python3 test_haznics.py
```

### Running demos, e.g. Poisson equation:
```
python3 poisson.py
```

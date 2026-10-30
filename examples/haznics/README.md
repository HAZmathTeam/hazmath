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
### Making `import haznics` work
The Python module is `swig_files/haznics.py` together with `swig_files/_haznics.so`,
built by `make install` with `haznics=yes`. Python only needs to find the
`swig_files` directory; `_haznics.so` finds `lib/libhazmath.so` on its own
(its runpath is set at build time), so `LD_LIBRARY_PATH` is not needed.
Choose one of the following.

- **Per shell session.** Run in bash, from any directory:
```
source /path/to/hazmath/examples/haznics/setup.rc
```
This adds `swig_files` to `PYTHONPATH` (once, even if sourced again) and runs
`dijitso clean` when legacy FEniCS is installed.

- **Every new shell.** Add the same `source` line to `~/.bashrc`.

- **Without environment variables.** Put a `.pth` file in the user
`site-packages` of the Python that will run the demos:
```
SITE=$(/usr/bin/python3 -m site --user-site)
mkdir -p "$SITE"
echo /path/to/hazmath/swig_files > "$SITE/haznics.pth"
```
That Python then finds `haznics` in every shell, script and IDE. Delete
`haznics.pth` to undo.

Check with
```
python3 -c "import haznics; print(haznics.__file__)"
```

Use the Python that matches the installed FEniCS/dolfinx and petsc4py (for the
Debian/Ubuntu packages, `/usr/bin/python3`); see `README_porting2dolfinx.md`.

### Running tests
```
cd tests/
python3 test_haznics.py
```

### Running demos, e.g. Poisson equation
```
python3 poisson.py
```

# haznicsx

HAZmath preconditioners for [DOLFINx](https://github.com/FEniCS/dolfinx) and
[PETSc](https://petsc.org), the successor of the SWIG-based `haznics` module
for legacy FEniCS.

This first version provides the core bindings:

- `CSRMatrix`: a HAZmath sparse matrix, built from CSR arrays
  (`csr_from_arrays`), a SciPy sparse matrix (`csr_from_scipy`) or an assembled
  sequential PETSc matrix (`csr_from_petsc`). It shares the arrays that need
  no conversion, so changing them afterwards changes the matrix.
- `AMGParam`: the HAZmath AMG parameters, with the HAZmath defaults.
- `amg(A, param)` and `famg(A, M, param)`: AMG and fractional AMG
  preconditioners. `Precond.apply(r, z)` computes `z = B r`.

PETSc `PC` classes and the auxiliary space (HX), rational approximation and
metric AMG preconditioners come in later versions.

## Installation

Both ways below run from the root of the HAZmath repository, which has the
`pyproject.toml`.

**In one go.** libhazmath is built and linked statically into the extension.
This needs a C/C++ compiler and `gfortran` (on Debian/Ubuntu:
`apt install build-essential gfortran`):

```bash
python3 -m pip install .[test]
```

HAZmath options are CMake definitions. For example, with SuiteSparse and
LAPACK (`apt install libsuitesparse-dev liblapack-dev`):

```bash
python3 -m pip install . -C cmake.define.USE_SUITESPARSE=ON -C cmake.define.USE_LAPACK=ON
```

**Against an installed libhazmath**, with its CMake package. Only a C/C++
compiler is needed for the second step:

```bash
cmake -S . -B build -DHAZMATH_RELOCATABLE_INSTALL=ON -DSHARED=ON \
      -DCMAKE_INSTALL_PREFIX=$HOME/.local/hazmath
cmake --build build --parallel
cmake --install build
CMAKE_PREFIX_PATH=$HOME/.local/hazmath python3 -m pip install .[test] \
    -C cmake.define.HAZNICSX_EXISTING_HAZMATH=ON
```

Either way, `src/haznics/helper.c` (the `create_precond_*` functions) is
compiled into the extension.

## Example

```python
import numpy as np
import scipy.sparse as sp

import haznicsx

n = 100
A = sp.diags([-1, 2, -1], [-1, 0, 1], shape=(n, n), format="csr")
B = haznicsx.amg(haznicsx.csr_from_scipy(A), haznicsx.AMGParam())

r = np.ones(n)
z = np.zeros(n)
B.apply(r, z)
```

## Limitations

- HAZmath uses 32-bit indices and real double precision values. Matrices
  from PETSc builds with complex or single precision scalars are rejected
  with a `TypeError`; indices from 64-bit PETSc builds are converted, and an
  `OverflowError` is raised if they do not fit.
- HAZmath is serial: `csr_from_petsc` accepts matrices on one process only.

## Tests

```bash
python3 -m pytest python/tests
```

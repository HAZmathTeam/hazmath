# haznicsx

HAZmath preconditioners for [DOLFINx](https://github.com/FEniCS/dolfinx) and
[PETSc](https://petsc.org), the successor of the SWIG-based `haznics` module
for legacy FEniCS.

It provides:

- `HazmathAMG` and `HazmathFAMG`: AMG and fractional AMG as petsc4py
  preconditioners (`PCPYTHON` contexts), configurable through the PETSc
  options database.
- `HazmathHXCurl` and `HazmathHXDiv`: HX (auxiliary space) preconditioners
  for H(curl) and H(div), in 2D and 3D. `hx_curl_operators(V)` and
  `hx_div_operators(V)` build their auxiliary matrices for a lowest order
  Nedelec or Raviart-Thomas space in DOLFINx.
- `CSRMatrix`: a HAZmath sparse matrix, built from CSR arrays
  (`csr_from_arrays`), a SciPy sparse matrix (`csr_from_scipy`), an assembled
  sequential PETSc matrix (`csr_from_petsc`) or a DOLFINx matrix
  (`csr_from_dolfinx`). It shares the arrays that need no conversion, so
  changing them afterwards changes the matrix.
- `AMGParam`: the HAZmath AMG parameters, with the HAZmath defaults.
- `amg`, `famg`, `hx_curl` and `hx_div`: the preconditioners without PETSc.
  `Precond.apply(r, z)` computes `z = B r`.

The rational approximation and metric AMG preconditioners come in later
versions.

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

HAZmath's index type `INT` is 32-bit by default. For problems with more than
about 2.1e9 indices, run `make config int64=yes` first; it writes
`include/haz_int.h`, which the build then uses. `haznicsx.index_dtype` tells
which width the package was built with.

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

## Examples

With DOLFINx, the preconditioner is chosen through the PETSc options:

```python
problem = dolfinx.fem.petsc.LinearProblem(
    a,
    L,
    bcs=[bc],
    petsc_options_prefix="poisson_",
    petsc_options={
        "ksp_type": "cg",
        "pc_type": "python",
        "pc_python_type": "haznicsx.HazmathAMG",
        "pc_hazmath_cycle_type": "w",  # any AMGParam field: pc_hazmath_<field>
    },
)
uh = problem.solve()
```

In parallel, use it on each process's block:
`"pc_type": "bjacobi", "sub_pc_type": "python", "sub_pc_python_type": "haznicsx.HazmathAMG"`.

The HX preconditioners need their auxiliary matrices, so they are attached
in code, e.g. for a Nedelec space `V`:

```python
ctx = haznicsx.HazmathHXCurl(*haznicsx.hx_curl_operators(V))
problem.solver.getPC().setType("python")
problem.solver.getPC().setPythonContext(ctx)
```

The additive variants (the default) are symmetric and work with CG; the
multiplicative ones (`multiplicative=True`) need e.g. GMRES.

With petsc4py directly:

```python
pc = ksp.getPC()
pc.setType("python")
pc.setPythonContext(haznicsx.HazmathAMG(haznicsx.AMGParam()))
```

Without PETSc:

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

- HAZmath uses real double precision values. Matrices from PETSc or DOLFINx
  with complex or single precision values are rejected with a `TypeError`.
  Indices of another width than `haznicsx.index_dtype` are converted, and an
  `OverflowError` is raised if they do not fit.
- HAZmath is serial: the converters and preconditioners accept matrices on
  one process only. In parallel, use `HazmathAMG` as the sub-preconditioner of
  block Jacobi or additive Schwarz; the HX preconditioners are serial only.
- A 64-bit `INT` needs the build in one go: the CMake install of libhazmath
  does not install `include/haz_int.h`, so code built against the installed
  headers assumes a 32-bit `INT`.
- PETSc objects that hold a haznicsx preconditioner should be destroyed
  before Python exits (as `LinearProblem` does); otherwise nanobind reports
  the haznicsx objects they still hold as leaked.

## Tests

```bash
python3 -m pytest python/tests
```

Tests that need PETSc are in files with `petsc` in the name, so they can be
deselected without PETSc: `python3 -m pytest python/tests --ignore-glob="*petsc*"`.

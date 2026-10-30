# Porting the HAZniCS Poisson demo to dolfinx

`poisson_dolfinx.py` is a dolfinx (0.11) version of `poisson.py`. It solves the
same problem, -div(grad(u)) + u = sin(pi*x) on a 16x16x16 unit cube with P1
elements and natural boundary conditions, with the same AMG parameters
(`max_levels=10`, `print_level=10`, `coarse_solver=32`).

## Running

```
PYTHONPATH=../../swig_files /usr/bin/python3 poisson_dolfinx.py
```

Use the Python that matches the installed dolfinx/petsc4py. The Debian/Ubuntu
packages of dolfinx, petsc4py and PETSc are built against the system numpy
(1.26.4). A Python environment that puts numpy 2.x first on the path stops at
the import of petsc4py with

```
ValueError: numpy.dtype size changed, may indicate binary incompatibility.
```

This is also why `poisson.py` (legacy dolfin + cbc.block, which imports the same
petsc4py) does not run in such an environment. Nothing in hazmath or haznics
needs numpy 2.x.

Result (serial): CG with hazmath UA AMG converges in 17 iterations to relative
tolerance 1e-10; the max difference from the PETSc sparse LU solution is
3.5e-10.

## What the dolfinx version does

cbc.block (`block.algebraic.hazmath.AMG`, `block.iterative.ConjGrad`) is not used.
Instead:

1. `A` and `b` are assembled with `dolfinx.fem.petsc.assemble_matrix` /
   `assemble_vector`.
2. The CSR arrays from `A.getValuesCSR()` are passed to `haznics.create_matrix`
   (values as float64, indices as int32, since hazmath `INT` is `int`).
3. The AMG is built with `haznics.AMG_param`, `param_amg_init`,
   `param_amg_set_dict` and `create_precond_amg`.
4. The AMG V-cycle is a PETSc `python` PC whose `apply` calls
   `haznics.apply_precond(r, z, prec)`; the outer solver is PETSc `cg`.
5. The comparison solution comes from PETSc `preonly` + `lu`, in place of
   dolfin's `solve(A, u, b)`.

The script is serial only: hazmath receives the whole matrix.

## Changes in the library needed for this

- `haz_shutils/headmk.sh`: source files named `__z*.c` are skipped when
  `include/hazmath.h` is generated. `src/fem/__zcommented_out_basis.c` and
  `src/utilities/__zcommented_out_io.c` contain only `#if 0` code, but their
  prototypes went into `hazmath.h` (14 functions, e.g. `bubble_face_basis`,
  `rt_basis`, `ned_basis`, `hazw`). SWIG wrapped them, and `import haznics`
  failed with `undefined symbol: bubble_face_basis`.
- `src/amr/scomplex.c` (two places) and `src/solver/amg_setup_ua.c` (one place):
  the `calloc` size is `MAX(n,0)`, which removes the three
  `-Walloc-size-larger-than=` warnings. Nothing changes for positive sizes.

The library was rebuilt with

```
make distclean ; make config suitesparse=yes shared=yes haznics=yes ; make install
```

without errors or warnings.

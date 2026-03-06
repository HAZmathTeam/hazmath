# Classical Ruge-Stuben AMG + Incomplete Factorization Preconditioners

This example demonstrates classical AMG and incomplete factorization
preconditioners from the HAZmath library.

## Library functions (in `src/amg_classical/`)

### Classical RS-AMG (`rs_classical.c`, `build_hierarchy.c`, `amg_cycle.c`)

- **`rs_amg_build_hierarchy`** -- Build multilevel AMG hierarchy from an SPD matrix
- **`rs_amg_rebuild_values`** -- Rebuild hierarchy values for a spectrally equivalent matrix (reuses coarsening)
- **`rs_amg_free`** -- Free all AMG data
- **`rs_amg_vcycle_precond`** -- Symmetric V-cycle: backslash + fwdslash
- **`rs_amg_backslash`** / **`rs_amg_fwdslash`** -- One-sided V-cycles
- **`rs_amg_ichol_precond`** -- Combined AMG + ichol preconditioner:
  `x = S*g + ichol(g - A*S*g) + S^T*(g - A*...)`

RS-AMG components: `rs_strength`, `rs_coarsening`, `rs_standard_interpolation`.

### Incomplete factorizations (`ichol.c`)

| Function | Description |
|---|---|
| `ichol_compute(A, L)` | Incomplete Cholesky, ichol(0), keeps sparsity of lower triangle |
| `icholt_compute(A, tau, L)` | Incomplete Cholesky with threshold dropping and diagonal compensation |
| `ichol_solve(L, b, x)` | Triangular solve: L L' x = b |
| `ilu_compute(A, L, U)` | ILU(0), keeps sparsity pattern of A |
| `ilut_compute(A, tau, L, U)` | ILU with threshold dropping, allows fill-in |
| `ilu_solve(L, U, b, x)` | Triangular solve: L U x = b |

All factorizations sort column indices via double transpose before factoring.

### Smoother types (from `macro.h`)

- `SMOOTHER_JACOBI` (1) -- damped Jacobi
- `SMOOTHER_GS` (2) -- lexicographic Gauss-Seidel
- `SMOOTHER_L1DIAG` (10) -- L1-Jacobi: D_l1 = diag(|A|*1)
- `SMOOTHER_GS_CF` (14) -- CF-ordered Gauss-Seidel

### AMG_param field mapping

| Field | Meaning | Default |
|---|---|---|
| `strong_coupled` | strength threshold | 0.25 |
| `coarse_dof` | min coarsest-level size | 4096 |
| `max_levels` | max AMG levels | 20 |
| `smoother` | smoother type | `SMOOTHER_GS_CF` |
| `presmooth_iter` | pre-smoothing sweeps | 2 |
| `postsmooth_iter` | post-smoothing sweeps | 2 |
| `relaxation` | Jacobi damping factor | 0.5 |

## Example programs

### test_amg

Tests AMG V-cycle and AMG+ichol as PCG preconditioners.

```
Usage: ./test_amg.ex <solve_matrix> [amg_matrix] [threshold] [nu] [min_size]
```

- Single matrix: builds AMG and solves with the same matrix
- Two matrices: builds AMG from `amg_matrix`, solves with `solve_matrix`

### test_ilu_ichol

Tests ichol(0), icholt, ILU(0), and ILUT as PCG preconditioners.

```
Usage: ./test_ilu_ichol.ex <matrix_file> [tau]
```

### amg_wrapper.c

Opaque-handle FFI wrapper for Julia/Python. Functions:
`amg_solver_create`, `amg_solver_solve`, `amg_solver_update`,
`amg_solver_set_matrix`, `amg_solver_free`, etc.

## Building

```bash
# Build HAZmath library first
make -C ../../ distclean
make -C ../../ config suitesparse=yes shared=yes
make -C ../../ install

# Build examples
make                              # builds test_amg.ex
make SRCFILE=test_ilu_ichol       # builds test_ilu_ichol.ex
make clean                        # removes all .ex and .o files
```

## Input format

Matrices are read in COO (coordinate) format via `dcoo_read_dcsr`:
```
nrow ncol nnz
i1 j1 val1
i2 j2 val2
...
```
Indices can be 0-based or 1-based (auto-detected). Matrices are
symmetrized as A = (A + A') / 2 before use.

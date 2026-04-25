# Classical RS-AMG and UA-AMG Solvers

This example demonstrates AMG preconditioners from the HAZmath library:
classical Ruge-Stuben (RS) AMG and unsmoothed aggregation (UA) AMG,
with incomplete factorization preconditioners.

## Library functions (in `src/amg_classical/`)

### Classical RS-AMG (`rs_classical.c`, `build_hierarchy.c`, `amg_cycle.c`)

- **`rs_amg_build_hierarchy`** -- Build multilevel AMG hierarchy from an SPD matrix
- **`rs_amg_rebuild_values`** -- Rebuild hierarchy values for a spectrally equivalent matrix (reuses coarsening)
- **`rs_amg_free`** -- Free all AMG data
- **`rs_amg_vcycle_precond`** -- Symmetric V/W-cycle: backslash + fwdslash
- **`rs_amg_backslash`** / **`rs_amg_fwdslash`** -- One-sided cycles
- **`rs_amg_ichol_precond`** -- Combined AMG + ichol preconditioner:
  `x = S*g + ichol(g - A*S*g) + S^T*(g - A*...)`

RS-AMG components: `rs_strength`, `rs_coarsening`, `rs_standard_interpolation`.

Coarsest-level solve uses UMFPACK direct solver (factorized once during hierarchy build).

Cycle type is controlled by `param.cycle_type`: 1 = V-cycle, 2 = W-cycle.

### UA-AMG (HAZmath built-in)

Uses HAZmath's `amg_setup_ua` for hierarchy construction (unsmoothed aggregation)
and `mgcycle` for the solve cycle. Default is W-cycle with Gauss-Seidel smoother.

### Incomplete factorizations (`ichol.c`)

| Function | Description |
|---|---|
| `ichol_compute(A, L)` | Incomplete Cholesky, ichol(0), keeps sparsity of lower triangle |
| `icholt_compute(A, tau, L)` | Incomplete Cholesky with threshold dropping and diagonal compensation |
| `ichol_solve(L, b, x)` | Triangular solve: L L' x = b |
| `ilu_compute(A, L, U)` | ILU(0), keeps sparsity pattern of A |
| `ilut_compute(A, tau, L, U)` | ILU with threshold dropping, allows fill-in |
| `ilu_solve(L, U, b, x)` | Triangular solve: L U x = b |

### Smoother types (from `macro.h`)

- `SMOOTHER_JACOBI` (1) -- damped Jacobi
- `SMOOTHER_GS` (2) -- lexicographic Gauss-Seidel
- `SMOOTHER_L1DIAG` (10) -- L1-Jacobi: D_l1 = diag(|A|*1)
- `SMOOTHER_GS_CF` (14) -- CF-ordered Gauss-Seidel

### RS-AMG parameter mapping

| Field | Meaning | Default |
|---|---|---|
| `strong_coupled` | strength threshold | 0.25 |
| `coarse_dof` | min coarsest-level size | 4096 |
| `max_levels` | max AMG levels | 20 |
| `smoother` | smoother type | `SMOOTHER_GS_CF` |
| `presmooth_iter` | pre-smoothing sweeps | 2 |
| `postsmooth_iter` | post-smoothing sweeps | 2 |
| `relaxation` | Jacobi damping factor | 0.5 |
| `cycle_type` | 1 = V-cycle, 2 = W-cycle | 1 |

## Example programs

### test_amg

Tests AMG as a PCG preconditioner. Supports two modes:

- **RS-AMG** (default): classical Ruge-Stuben with V or W-cycle,
  UMFPACK direct solve on the coarsest level.
- **UA-AMG**: unsmoothed aggregation with W-cycle (set `USE_UA=1`).

```
Usage: ./test_amg.ex [--cycle=V|W|UA] <solve_matrix> [amg_matrix] [threshold] [nu] [min_size]

  --cycle=V   RS AMG, V-cycle (default)
  --cycle=W   RS AMG, W-cycle
  --cycle=UA  UA AMG, W-cycle

Examples:
  ./test_amg.ex matrix.txt                        # RS V-cycle
  ./test_amg.ex --cycle=W matrix.txt              # RS W-cycle
  ./test_amg.ex --cycle=UA matrix.txt             # UA W-cycle
  ./test_amg.ex solve.txt amg.txt 0.25 2 4096     # two-matrix RS V-cycle
```

- Single matrix: builds AMG and solves with the same matrix
- Two matrices: builds AMG from `amg_matrix`, solves with `solve_matrix`

Sample results on a 26M DOF problem:

| Method | Levels | Setup | Iter | Solve | Total | RAM |
|---|---|---|---|---|---|---|
| RS V-cycle | 7 | 61s | 9 | 73s | 134s | 18.6 GB |
| RS W-cycle | 7 | 63s | 2 | 69s | 132s | 18.7 GB |
| UA W-cycle | 11 | 79s | 18 | 304s | 383s | 10.5 GB |

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
# Build HAZmath library first (suitesparse required for UMFPACK coarsest-level solve)
make -C ../../ distclean
make -C ../../ config suitesparse=yes lapack=yes
make -C ../../ install

# Build examples
make                              # builds test_amg.ex and test_ilu_ichol.ex
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

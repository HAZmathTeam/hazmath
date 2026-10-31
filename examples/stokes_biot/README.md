# stokes_biot: Bernardi-Raugel / P0 Stokes

Solves

    -div(2 mu eps(u)) + grad p = f,   div u = g   in (0,1)^d,  d = 2, 3,
    u = 0 on Gamma_D,   (2 mu eps(u) - p I) n = g_N on Gamma_N,

with the Bernardi-Raugel element (P1 velocity + normal face bubbles) and P0
pressure, on a sequence of uniform meshes, and prints the errors against a
manufactured solution and the convergence rates.  The Biot problem will be
added here later.

## Build

The HAZmath library must be built with SuiteSparse and LAPACK
(`make config suitesparse=yes lapack=yes shared=yes; make install` in the
HAZmath root).  Then

    make            # -> stokes_br.ex

## Run

    ./stokes_br.ex [mu=<value>] [test=0|1|2] [rhs=pr|std] [condensed|full]
                   [solver=gmres|umfpack|compare] [amg=rs|ua] [cycle=amli|w|v]

The options are read from `input.dat` (AMG), then from `stokes_br.dat`
(options specific to this example), then from the command line, each
overriding the previous.  With the files as shipped and no arguments: mu = 1,
test 0, pressure-robust right-hand side, condensed bubbles, GMRES with the
RS-AMG block preconditioner.

| option | file key | values (default first) | meaning |
|---|---|---|---|
| `mu=` | `viscosity` (stokes_br.dat) | `1` | viscosity |
| `test=` | `test_problem` (stokes_br.dat) | `0`, `1`, `2` | 0: div u = 0, Dirichlet on the whole boundary; 1: div u != 0, Dirichlet on the whole boundary; 2: div u != 0, Dirichlet on x_1 = 0, traction on the rest |
| `rhs=` | `rhs` (stokes_br.dat) | `pr`, `std` | `pr`: source and normal traction tested with the BDM1 interpolant Pi_h v (pressure robust); `std`: tested with v |
| `condensed`, `full` | `system` (stokes_br.dat) | `condensed`, `full` | `condensed`: static condensation of the bubbles (bubble-bubble block replaced by its diagonal; only the diagonal is assembled); `full`: the full BR/P0 system |
| `solver=` | `solver` (stokes_br.dat) | `gmres`, `umfpack`, `compare` | `gmres`: right-preconditioned GMRES, tol 1e-10; if it fails with RS-AMG, GMRES with UA-AMG (AMLI-cycle) is tried, then UMFPACK. `compare`: GMRES and UMFPACK, last column is their relative difference |
| `amg=` | `AMG_type` (input.dat) | `rs` (`C`), `ua` (`UA`) | AMG in the block lower triangular preconditioner: classical Ruge-Stuben (V-cycle) or unsmoothed aggregation |
| `cycle=` | `AMG_cycle_type` (input.dat) | `amli` (`A`), `w` (`W`), `v` (`V`) | cycle for UA only (RS always uses the V-cycle) |

Mesh and quadrature parameters are in `input.dat`:

| parameter | meaning |
|---|---|
| `spatial_dim` | 2 or 3 |
| `refinement_levels` | the program solves on levels 1, ..., refinement_levels (h = 2^-level) |
| `refinement_type` | > 10: uniform refinement |
| `nquad` | quadrature points per direction, on the elements and on the traction faces (test 2); with `rhs=pr` and small `mu` use `nquad = 5`. In 3D the element rule (`tetquad_`) is wrong for `nquad > 5` |

Examples:

    ./stokes_br.ex                          # defaults
    ./stokes_br.ex mu=1e-6 test=2           # traction BC, small viscosity
    ./stokes_br.ex solver=compare amg=ua    # GMRES (UA, AMLI) against UMFPACK

The output is one line per level: number of dofs (full BR system), ||u-u_h||,
|u-u_h|_1, ||p-p_h|| (modulo constants when Gamma_N is empty),
max_T |int_T div(u-u_h)|, GMRES iterations and CPU times.  Long runs: redirect
to a log, e.g. `./stokes_br.ex mu=1e-6 2>/dev/null | tee output/run.log`.

## Files

| file | content |
|---|---|
| `stokes_br.c` | driver: mesh, spaces, assembly, condensation, solve, errors |
| `stokes_br_data.h` | test problems: exact solutions, f, g, traction, Gamma_N |
| `stokes_br_system.h` | local matrix and right-hand side, bubble condensation and recovery, errors |
| `stokes_br_precond.h` | block lower triangular AMG preconditioner, GMRES and UMFPACK solves |
| `doc/stokes_br_pr.tex` | note on the pressure-robust right-hand sides |
| `input.dat` | mesh, quadrature and AMG parameters (read by the HAZmath `param_input`) |
| `stokes_br.dat` | options specific to this example (viscosity, test problem, rhs, condensed/full, solver) |
| `regression.sh` | runs 8 cases on levels 1-4 and compares with `reference/`; `--update` regenerates the references |
| `reference/` | reference outputs (timing columns removed) |

## Regression tests

`./regression.sh` prints, for each case, PASS/FAIL, the wall time and the
peak memory (maximum resident set size) of the run, levels 1-4, nquad = 5.
On an Intel Core Ultra X7 358H:

| case | options | wall (s) | memory (MB) |
|---|---|---|---|
| `2d_default` | | 0.07 | 11.0 |
| `2d_test1` | `test=1` | 0.05 | 10.8 |
| `2d_test2_mu1e-6` | `test=2 mu=1e-6` | 0.06 | 11.2 |
| `2d_full_umfpack` | `full solver=umfpack` | 0.08 | 19.3 |
| `2d_ua_amli` | `amg=ua` | 0.06 | 10.0 |
| `3d_default` | | 4.24 | 237.4 |
| `3d_test2_mu1e-6` | `test=2 mu=1e-6` | 7.09 | 238.9 |
| `3d_ua_amli` | `amg=ua` | 4.81 | 223.4 |

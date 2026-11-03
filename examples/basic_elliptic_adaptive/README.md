# Adaptive finite elements for grad, curl and div problems (HAZmath)

Adaptive finite elements on the L-shaped domain (2D) and the Fichera domain
(3D). As in `examples/basic_elliptic`, `FE_type` in `input.dat` selects the
element:

| `FE_type` | element | problem | estimator |
|---|---|---|---|
| 1 | P1 | `-Delta u = f` | by preconditioning: local problems on vertex patches (`eta`) and the explicit residual estimator (`zeta`), [1] |
| 20 | lowest order Nedelec | `(eps curl u, curl v) + (kappa u, v) = (f, v)`, `u x n = 0` | robust residual estimator, [2, Thm 5.8, (5.19)] |
| 30 | lowest order Raviart-Thomas | `(eps div u, div v) + (kappa u, v) = (f, v)`, `u . n = 0` | robust residual estimator, [2, Thm 5.8, (5.20)] |

[1] Y. Li and L. Zikatanov, *A posteriori error estimates of finite element
methods by preconditioning*, Comput. Math. Appl. 91 (2021) 192-201,
doi:10.1016/j.camwa.2020.08.001.

[2] Y. Li and L. Zikatanov, *Nodal auxiliary a posteriori error estimates*,
Math. Comp. 96 (2027) 47-75, doi:10.1090/mcom/4141.

The adaptive loop is SOLVE -> ESTIMATE -> MARK -> REFINE:

- SOLVE: CG with AMG (P1), or CG with the HX preconditioner for Nedelec
  (curl) and Raviart-Thomas (div, 3D) with 10000 dofs or more; the sparse
  direct solver for smaller Nedelec/Raviart-Thomas systems and for
  Raviart-Thomas in 2D;
- ESTIMATE: element indicators whose sum is the estimator;
- MARK: Doerfler marking with parameter `ae_theta` on the squared indicators
  (`ae_theta{1.0}` marks every element, i.e. uniform bisection);
- REFINE: newest vertex bisection (Diening-Gehring-Storn).

## Build

```sh
cd HAZMATH_DIR
make distclean; make config suitesparse=yes lapack=yes shared=yes; make install
cd examples/basic_elliptic_adaptive
make clean; make
```

## Run

Two files control a run:

- `input.dat` (read from the current directory): `FE_type` and the solver
  parameters;
- the grid file given on the command line (`input/2d_L.input` or
  `input/3d_fichera.input`): the domain, `ae_problem`, `ae_eps`, `ae_kappa`,
  `ae_theta` and the other adaptivity parameters (table below).

```sh
./basic_elliptic_adaptive.ex input/2d_L.input        # L-shaped domain (2D)
./basic_elliptic_adaptive.ex input/3d_fichera.input  # Fichera domain (3D)
```

### Poisson (P1)

In `input.dat`:

    FE_type = 1

In the grid file: `ae_problem{1}` (L-shape only, singular exact solution),
`ae_problem{2}` (`f = 1`) or `ae_problem{3}` (smooth exact solution).
The files as shipped run P1: `2d_L.input` has `ae_problem{1}`,
`3d_fichera.input` has `ae_problem{2}`.

### Curl-curl (Nedelec)

In `input.dat`:

    FE_type = 20

In the grid file: `ae_problem{3}` (smooth exact solution, the effectivity
index `eta/err` is printed) or `ae_problem{2}` (`f = (1,...,1)`);
`ae_problem{1}` is P1 only, so change it in `2d_L.input`. Set the
coefficients with `ae_eps{...}` and `ae_kappa{...}` (default 1 and 1), e.g.

    ae_problem{3}
    ae_eps{0.1}
    ae_kappa{10.}

### Div-div (Raviart-Thomas)

In `input.dat`:

    FE_type = 30

The grid file keys are the same as for curl-curl (`ae_problem{3}` or `{2}`,
`ae_eps`, `ae_kappa`).

The convergence table is printed on the screen and written to `output/`
(see Output below). For long 3D runs `ae_max_dofs{N}` stops the loop after the
first mesh with more than `N` dofs.

## Problems

| `ae_problem` | `FE_type` 1 | `FE_type` 20, 30 |
|---|---|---|
| 1 | 2D only: `f = 0`, exact `u = r^{2/3} sin(2/3 (theta - pi/2))` | not available |
| 2 | `f = 1`, `u = 0` on the boundary | `f = (1,...,1)` |
| 3 | exact `u = prod_i sin(pi x_i)` | exact: curl 2D `u = (cos(pi x) sin(pi y), sin(pi x) cos(pi y))`, curl 3D `u = (0, 0, sin(pi x) sin(pi y))`, div `u = (sin(pi x_1), ..., sin(pi x_d))` |

Problem 2 has no closed-form solution; its error is computed from
`||u - u_h||^2 = (f,u) - (f,u_h)` when `ae_reference_energy` gives `(f,u)`.

For P1 with constant `f` (problems 1 and 2) everything is integrated exactly.
Quadrature (`nquad` in `input.dat`) is used for the right-hand side and the
error with a non-constant `f`, and for the Nedelec and Raviart-Thomas
right-hand sides. The Nedelec and Raviart-Thomas matrices are computed exactly.

## Parameters in the grid file

| key | meaning |
|---|---|
| `num_refinements{n}` | maximal number of stages of the adaptive loop |
| `ae_problem{p}` | table above |
| `ae_theta{t}` | Doerfler parameter; 1 = uniform |
| `ae_estimator{e}` | P1 only: marking with 0 = `eta` (local problems), 1 = `zeta` (residual) |
| `ae_reference_energy{E}` | problem 2: `(f,u)`; `<= 0` = unknown |
| `ae_max_dofs{N}` | stop after the first mesh with more than `N` dofs |
| `ae_eps{e}`, `ae_kappa{k}` | `FE_type` 20, 30: the coefficients |
| `ae_dump_meshes{0/1}` | 1: write `output/mesh_stage_XXX.vtu` at every stage |

## Output

- P1: `output/convergence_problem<p>.txt`: one line per stage with the number
  of dofs and elements, the energy error, `eta`, `zeta`, the effectivity
  indices `eta/err` and `zeta/err`, `(f,u_h)`, the number of marked elements
  and the oscillation; `output/sol_problem<p>.vtu` and
  `output/mesh_problem<p>.msh` at the last stage;
- `FE_type` 20, 30: `output/convergence_fe<FE_type>_problem<p>.txt`: one line
  per stage with the number of dofs and elements, the solver iterations
  (0 = direct), the energy error, `eta`, `eta/err`, `(f,u_h)` and the number
  of marked elements.

The derivation of the P1 estimators is in `doc/basic_elliptic_adaptive.tex`.

## Animated GIF of the refinement (optional)

With `ae_dump_meshes{1}` in the grid file, every stage writes
`output/mesh_stage_XXX.vtu`. With ParaView (`pvpython`) and ImageMagick:

```sh
pvpython --force-offscreen-rendering render_stages.py output frames "Fichera, P1 Laplace, f = 1"
cd frames
for f in left_*.png; do i=${f#left_}; convert $f right_$i +append +repage -resize 1300x +repage both_$i; done
convert -delay 70 both_*.png -loop 0 -layers Optimize ../refinement.gif
```

`render_stages.py` draws two panels per stage (the whole domain and a zoom near
the re-entrant corner of the Fichera domain). Smaller files:
`-resize 900x +repage -fuzz 3% -colors 32` before `-loop 0`. The `.vtu` files
are large (about 250 MB for 20 stages up to 1M elements); delete them after
making the GIF.

Created by Ludmil Zikatanov, 2026-11-03, with the help of Claude (Anthropic).

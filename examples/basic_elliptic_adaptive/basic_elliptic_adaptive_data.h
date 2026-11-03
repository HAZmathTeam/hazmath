/*! \file examples/basic_elliptic_adaptive/basic_elliptic_adaptive_data.h
 *
 *  Created by Ludmil Zikatanov 20261103.
 *  Copyright 2026__HAZMATH__. All rights reserved.
 *  Created with the help of Claude (Anthropic).
 *
 * \brief Data for the adaptive P1 Laplace example on the L-shaped domain
 *        Omega = (-1,1)^2 \ [0,1]^2 (re-entrant corner at the origin) and
 *        on the Fichera domain (-1,1)^3 \ [0,1]x[-1,0]x[0,1].
 *
 *   problem 1:  (2D only) f = 0, u = g on the boundary, exact solution
 *               u = r^{2/3} sin(2/3 (theta - pi/2)),  theta in [pi/2, 2 pi].
 *   problem 2:  f = 1, u = 0 on the boundary (no closed-form solution).
 *   problem 3:  u = prod_i sin(pi x_i), f = d pi^2 u (2D and 3D).
 *
 * \note In problems 1 and 2 f is constant, so the P1 load vector and the
 *       estimators are computed exactly, without quadrature (see doc/).
 */

/* angle in [pi/2, 2 pi] on the closure of the L-shaped domain */
static REAL lshape_theta(REAL* x) {
  REAL t = atan2(x[1], x[0]);
  if (t <= 0.0) t += 2.0 * M_PI;
  return t;
}

/* problem 1: exact solution */
void exactsol_lshape(REAL* val, REAL* x, REAL time, void* param) {
  REAL r = sqrt(x[0] * x[0] + x[1] * x[1]);
  if (r < 1e-15) { *val = 0.0; return; }
  *val = pow(r, 2.0 / 3.0) * sin(2.0 / 3.0 * (lshape_theta(x) - 0.5 * M_PI));
}

/* problem 1: gradient of the exact solution,
   grad u = (2/3) r^{-1/3} ( -sin(theta/3 + pi/3), cos(theta/3 + pi/3) ) */
void D_exactsol_lshape(REAL* val, REAL* x, REAL time, void* param) {
  REAL r = sqrt(x[0] * x[0] + x[1] * x[1]);
  if (r < 1e-15) { val[0] = 0.0; val[1] = 0.0; return; }
  REAL a = lshape_theta(x) / 3.0 + M_PI / 3.0;
  REAL c = (2.0 / 3.0) * pow(r, -1.0 / 3.0);
  val[0] = -c * sin(a);
  val[1] = c * cos(a);
}

/* Dirichlet data */
void bc_lshape_p1(REAL* val, REAL* x, REAL time, void* param) {
  exactsol_lshape(val, x, time, param);
}
void bc_lshape_p2(REAL* val, REAL* x, REAL time, void* param) {
  *val = 0.0;
}

/* constant right-hand side for problems 1 and 2 */
static REAL rhs_const_lshape(INT problem) {
  return (problem == 2) ? 1.0 : 0.0;
}

/* problem 3: u = prod_i sin(pi x_i), which vanishes on the boundary of the
   L-shape and of the Fichera domain (all faces lie on x_i in {-1,0,1});
   f = -Delta u = d pi^2 u. The dimension is set in main(). */
static INT g_dim_sines = 2;
void exactsol_sines(REAL* val, REAL* x, REAL time, void* param) {
  REAL s = 1.0;
  for (INT i = 0; i < g_dim_sines; i++) s *= sin(M_PI * x[i]);
  *val = s;
}
void D_exactsol_sines(REAL* val, REAL* x, REAL time, void* param) {
  for (INT j = 0; j < g_dim_sines; j++) {
    REAL s = M_PI;
    for (INT i = 0; i < g_dim_sines; i++) s *= (i == j) ? cos(M_PI * x[i]) : sin(M_PI * x[i]);
    val[j] = s;
  }
}
void rhs_sines(REAL* val, REAL* x, REAL time, void* param) {
  exactsol_sines(val, x, time, param);
  *val *= (REAL)g_dim_sines * M_PI * M_PI;
}

/*********************************************************************/
/* FE_type 20 (Ned0) and 30 (RT0):
 *   curl: (eps curl u, curl v) + (kappa u, v) = (f, v),  u x n = 0,
 *   div:  (eps div u, div v)   + (kappa u, v) = (f, v),  u . n = 0.
 *
 *   problem 2: f = (1,...,1), no closed-form solution.
 *   problem 3: smooth exact solution; it satisfies the boundary condition
 *              on every domain whose faces lie on x_i in {-1, 0, 1}:
 *     curl, 2D: u = (cos(pi x) sin(pi y), sin(pi x) cos(pi y)),  f = kappa u,
 *     curl, 3D: u = (0, 0, sin(pi x) sin(pi y)),  f = (2 pi^2 eps + kappa) u,
 *     div:      u = (sin(pi x_1), ..., sin(pi x_d)),  f = (pi^2 eps + kappa) u.
 */
static void hd_exact_u(INT fet, INT dim, const REAL* x, REAL* u) {
  INT j;
  if (fet == 20 && dim == 2) {
    u[0] = cos(M_PI * x[0]) * sin(M_PI * x[1]);
    u[1] = sin(M_PI * x[0]) * cos(M_PI * x[1]);
  } else if (fet == 20) {
    u[0] = 0.0; u[1] = 0.0;
    u[2] = sin(M_PI * x[0]) * sin(M_PI * x[1]);
  } else {
    for (j = 0; j < dim; j++) u[j] = sin(M_PI * x[j]);
  }
}
/* D u: curl u (1 component in 2D, 3 in 3D) or div u (1 component) */
static void hd_exact_du(INT fet, INT dim, const REAL* x, REAL* d) {
  if (fet == 20 && dim == 2) {
    d[0] = 0.0;
  } else if (fet == 20) {
    d[0] = M_PI * sin(M_PI * x[0]) * cos(M_PI * x[1]);
    d[1] = -M_PI * cos(M_PI * x[0]) * sin(M_PI * x[1]);
    d[2] = 0.0;
  } else {
    d[0] = 0.0;
    for (INT j = 0; j < dim; j++) d[0] += M_PI * cos(M_PI * x[j]);
  }
}
static void hd_rhs(INT fet, INT dim, INT problem, REAL eps, REAL kap, const REAL* x, REAL* f) {
  INT j;
  if (problem == 2) { for (j = 0; j < dim; j++) f[j] = 1.0; return; }
  hd_exact_u(fet, dim, x, f);
  REAL s;
  if (fet == 20) s = (dim == 2) ? kap : (2.0 * M_PI * M_PI * eps + kap);
  else s = M_PI * M_PI * eps + kap;
  for (j = 0; j < dim; j++) f[j] *= s;
}
/* the data part of the element residual R1:
   curl: -div f (1 component); div: curl f (1 component in 2D, 3 in 3D).
   (div of Ned0 and curl of RT0 functions vanish on each element.) */
static void hd_R1(INT fet, INT dim, INT problem, REAL eps, REAL kap, const REAL* x, REAL* r) {
  r[0] = 0.0; r[1] = 0.0; r[2] = 0.0;
  if (problem == 3 && fet == 20 && dim == 2)
    r[0] = 2.0 * M_PI * kap * sin(M_PI * x[0]) * sin(M_PI * x[1]);
}

/*! \file examples/stokes_biot/stokes_br_data.h
 *
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * \brief Problem data for the Bernardi-Raugel Stokes example:
 *
 *        -div(2 mu eps(u)) + grad(p) = f,   div(u) = g   in (0,1)^d,
 *        u = 0 on Gamma_D,   (2 mu eps(u) - p I) n = g_N on Gamma_N.
 *
 *        Test problems (stokes_test):
 *        0: a = t^2 (1-t)^2,
 *             2D: psi = a(x) a(y),        u = (d_y psi, -d_x psi),
 *             3D: phi = a(x) a(y) a(z),   u = curl(phi, phi, phi),
 *           div u = 0, Gamma_D = boundary.
 *        1: u_i = A_i prod_k sin(pi x_k), div u != 0, Gamma_D = boundary.
 *        2: u_i = A_i sin(pi x_1/2) prod_{k>1} sin(pi x_k + pi/4), div u != 0,
 *           Gamma_D = {x_1 = 0}, traction on the rest of the boundary.
 *        In all tests p = prod_k cos(pi x_k) and
 *          f = -mu (Lap(u) + grad(div u)) + grad(p),  g = div u,
 *          g_N = mu (grad u + grad u^T) n - p n.
 *
 * \note Block layout of the unknowns: [ bubble, u_1, ..., u_d, p ].
 *
 */

// Viscosity and test problem (set from the command line in stokes_br.c)
REAL stokes_mu = 1.0;
INT stokes_test = 0;

// a(t) = t^2 (1-t)^2 and its derivatives of order m = 0,1,2,3
static REAL a_der(REAL t, INT m)
{
  switch (m) {
  case 0:  return t * t * (1.0 - t) * (1.0 - t);
  case 1:  return 2.0 * t - 6.0 * t * t + 4.0 * t * t * t;
  case 2:  return 2.0 - 12.0 * t + 12.0 * t * t;
  case 3:  return -12.0 + 24.0 * t;
  default: return 0.0;
  }
}

// d^{m_0}_x d^{m_1}_y (d^{m_2}_z) of the potential prod_k a(x_k)
static REAL dpot(REAL* x, INT dim, INT m0, INT m1, INT m2)
{
  REAL val = a_der(x[0], m0) * a_der(x[1], m1);
  if (dim == 3) val *= a_der(x[2], m2);
  return val;
}

// Derivative of the potential: multi-index = e_i + e_j + e_k (i,j,k = -1 for none)
static REAL dpot_ijk(REAL* x, INT dim, INT i, INT j, INT k)
{
  INT m[3] = {0, 0, 0};
  if (i >= 0) m[i]++;
  if (j >= 0) m[j]++;
  if (k >= 0) m[k]++;
  return dpot(x, dim, m[0], m[1], m[2]);
}

// Tests 1, 2: d_j d_k u_i for u_i = A_i prod_l sin(w_l x_l + s_l)
static REAL u_der_trig(REAL* x, INT dim, INT i, INT j, INT k)
{
  const REAL A[3] = {1.0, -0.5, 0.75};
  INT l, m[3] = {0, 0, 0};
  if (j >= 0) m[j]++;
  if (k >= 0) m[k]++;
  REAL val = A[i];
  for (l = 0; l < dim; l++) {
    REAL w = M_PI, s = 0.0;
    if (stokes_test == 2) {
      if (l == 0) w = 0.5 * M_PI;
      else s = 0.25 * M_PI;
    }
    val *= pow(w, m[l]) * sin(w * x[l] + s + 0.5 * M_PI * m[l]);
  }
  return val;
}

// d_j d_k u_i (j,k = -1 for none).  Test 0: u = (d_y psi, -d_x psi) in 2D,
// u_i = d_{i+1} phi - d_{i+2} phi in 3D.
static REAL u_der(REAL* x, INT dim, INT i, INT j, INT k)
{
  if (stokes_test > 0) return u_der_trig(x, dim, i, j, k);
  if (dim == 2) {
    if (i == 0) return  dpot_ijk(x, dim, 1, j, k);
    else        return -dpot_ijk(x, dim, 0, j, k);
  }
  return dpot_ijk(x, dim, (i + 1) % 3, j, k) - dpot_ijk(x, dim, (i + 2) % 3, j, k);
}

// Exact velocity u
void exact_u(REAL* val, REAL* x, INT dim)
{
  INT i;
  for (i = 0; i < dim; i++) val[i] = u_der(x, dim, i, -1, -1);
}

// Exact velocity gradient: val[i*dim + j] = d_j u_i
void Dexact_u(REAL* val, REAL* x, INT dim)
{
  INT i, j;
  for (i = 0; i < dim; i++)
    for (j = 0; j < dim; j++) val[i * dim + j] = u_der(x, dim, i, j, -1);
}

// Exact divergence g = div u
REAL exact_div(REAL* x, INT dim)
{
  INT j;
  REAL g = 0.0;
  for (j = 0; j < dim; j++) g += u_der(x, dim, j, j, -1);
  return g;
}

// Exact pressure p = prod_k cos(pi x_k)
REAL exact_p(REAL* x, INT dim)
{
  INT k;
  REAL p = 1.0;
  for (k = 0; k < dim; k++) p *= cos(M_PI * x[k]);
  return p;
}

// Momentum source f = -div(2 mu eps(u)) + grad(p) = -mu (Lap(u) + grad(div u)) + grad(p)
void source_f(REAL* val, REAL* x, INT dim)
{
  INT i, j, k;
  for (i = 0; i < dim; i++) {
    val[i] = 0.0;
    for (j = 0; j < dim; j++)
      val[i] -= stokes_mu * (u_der(x, dim, i, j, j) + u_der(x, dim, j, i, j));
    REAL dp = -M_PI * sin(M_PI * x[i]);
    for (k = 0; k < dim; k++)
      if (k != i) dp *= cos(M_PI * x[k]);
    val[i] += dp;
  }
}

// Traction g_N = (2 mu eps(u) - p I) n
void traction(REAL* val, REAL* x, REAL* n, INT dim)
{
  INT i, j;
  REAL G[9];
  Dexact_u(G, x, dim);
  REAL p = exact_p(x, dim);
  for (i = 0; i < dim; i++) {
    val[i] = -p * n[i];
    for (j = 0; j < dim; j++) val[i] += stokes_mu * (G[i * dim + j] + G[j * dim + i]) * n[j];
  }
}

/* hazmath callbacks.  RHS layout follows FEM_Block_RHS_Local:
 * [ bubble (dim entries), u_1, ..., u_d, p ] -> 2*dim+1 entries. */
INT stokes_br_dim = 2;

void source_rhs(REAL* val, REAL* x, REAL time, void* param)
{
  INT k, dim = stokes_br_dim;
  REAL f[3];
  source_f(f, x, dim);
  for (k = 0; k < dim; k++) val[k] = f[k];        // bubble test functions
  for (k = 0; k < dim; k++) val[dim + k] = f[k];  // P1 test functions
  val[2 * dim] = -exact_div(x, dim);              // -(div u, q) = -(g, q)
}

// Homogeneous Dirichlet data for all unknowns
void bc_zero(REAL* val, REAL* x, REAL time, void* param)
{
  INT k;
  for (k = 0; k < 2 * stokes_br_dim + 1; k++) val[k] = 0.0;
}

// Is x on the traction boundary Gamma_N (test 2: boundary minus {x_1 = 0})?
INT on_gamma_N(REAL* x, INT dim)
{
  INT k;
  if (stokes_test != 2) return 0;
  if (x[0] < 1e-12) return 0;
  for (k = 0; k < dim; k++)
    if (x[k] < 1e-12 || x[k] > 1.0 - 1e-12) return 1;
  return 0;
}

/*
 * amg_cycle.c - AMG backslash/fwdslash cycles and V-cycle preconditioner
 *
 * Authors: HAZmath (https://hazmath.net)
 * Created with the help of Claude (Anthropic)
 *
 * Smoothers: L1-Jacobi (default), Gauss-Seidel, CF-ordered GS, damped Jacobi.
 * Uses HAZmath smoother functions where possible.
 * Coarsest level solved by ichol-PCG.
 */
#include "hazmath.h"

/* ichol preconditioner wrapper for coarsest-level solve */
static void ichol_precond_coarse(REAL* r, REAL* z, void* data) {
  const dCSRmat* L = (const dCSRmat*)data;
  ichol_solve(L, r, z);
}

/* Coarsest-level solve: direct (UMFPACK) or ichol-PCG */
static void coarse_solve(dCSRmat* A, const REAL* b, REAL* x,
                         const dCSRmat* L_ichol, INT n) {
  dvector bv = {n, (REAL*)b};
  dvector uv = {n, x};
  if (L_ichol && L_ichol->nnz > 0) {
    /* ichol-preconditioned PCG (for SPD matrices) */
    precond pc;
    pc.data = (void*)L_ichol;
    pc.fct  = ichol_precond_coarse;
    dcsr_pcg(A, &bv, &uv, &pc, 1e-12, n, STOP_REL_PRECRES, 0);
  } else {
    /* Direct solve via UMFPACK (works for non-symmetric) */
    directsolve_HAZ(A, &bv, &uv, 0);
  }
}

/* ======================================================================
 *  L1-Jacobi smoother (Baker, Falgout, Kolev, Yang 2011)
 *
 *    x = x + D_l1^{-1} (b - Ax)
 *    (D_l1)_ii = sum_j |a_ij|  (l1 row norm)
 * ====================================================================== */
static void l1_smooth(dCSRmat* A, const REAL* l1_diag, INT n,
                      const REAL* b, REAL* x, INT nu) {
  REAL* r = (REAL*)malloc(n * sizeof(REAL));
  for (INT sweep = 0; sweep < nu; sweep++) {
    /* r = b - A*x */
    array_cp(n, (REAL*)b, r); dcsr_aAxpy(-1.0, A, x, r);
    /* x = x + D_l1^{-1} * r */
    for (INT i = 0; i < n; i++)
      x[i] += r[i] / l1_diag[i];
  }
  free(r);
}

/* ======================================================================
 *  CF-ordered Gauss-Seidel
 * ====================================================================== */
static void cf_fwd_gs(dCSRmat* A, const INT* cf_order, INT n,
                      const REAL* b, REAL* x, INT nu) {
  for (INT sweep = 0; sweep < nu; sweep++) {
    for (INT ii = 0; ii < n; ii++) {
      INT i = cf_order[ii];
      REAL rhs = b[i], diag = 1.0;
      for (INT k = A->IA[i]; k < A->IA[i + 1]; k++) {
        if (A->JA[k] == i) diag = A->val[k];
        else               rhs -= A->val[k] * x[A->JA[k]];
      }
      x[i] = rhs / diag;
    }
  }
}

static void cf_bwd_gs(dCSRmat* A, const INT* cf_order, INT n,
                      const REAL* b, REAL* x, INT nu) {
  for (INT sweep = 0; sweep < nu; sweep++) {
    for (INT ii = n - 1; ii >= 0; ii--) {
      INT i = cf_order[ii];
      REAL rhs = b[i], diag = 1.0;
      for (INT k = A->IA[i]; k < A->IA[i + 1]; k++) {
        if (A->JA[k] == i) diag = A->val[k];
        else               rhs -= A->val[k] * x[A->JA[k]];
      }
      x[i] = rhs / diag;
    }
  }
}

/* ======================================================================
 *  Smoother dispatch: forward (pre-smooth) and backward (post-smooth)
 *
 *  Uses hazmath smoothers for GS and Jacobi, custom for L1 and CF-GS.
 * ====================================================================== */
static void smooth_fwd(AMG_data* mgl, AMG_param* param, INT lev,
                       const REAL* b, REAL* x, INT nu) {
  rs_level_aux* aux = RS_AUX(mgl, lev);
  INT n = mgl[lev].A.row;

  switch (param->smoother) {
    case SMOOTHER_GS: {
      dvector u_vec = {n, x};
      dvector b_vec = {n, (REAL*)b};
      smoother_dcsr_gs(&u_vec, 0, n-1, 1, &mgl[lev].A, &b_vec, nu);
      break;
    }
    case SMOOTHER_GS_CF:
      cf_fwd_gs(&mgl[lev].A, aux->cf_order, n, b, x, nu);
      break;
    case SMOOTHER_JACOBI: {
      dvector u_vec = {n, x};
      dvector b_vec = {n, (REAL*)b};
      smoother_dcsr_jacobi(&u_vec, 0, n-1, 1, &mgl[lev].A, &b_vec, nu);
      break;
    }
    case SMOOTHER_L1DIAG:
    default:
      l1_smooth(&mgl[lev].A, aux->l1_diag, n, b, x, nu);
      break;
  }
}

static void smooth_bwd(AMG_data* mgl, AMG_param* param, INT lev,
                       const REAL* b, REAL* x, INT nu) {
  rs_level_aux* aux = RS_AUX(mgl, lev);
  INT n = mgl[lev].A.row;

  switch (param->smoother) {
    case SMOOTHER_GS: {
      dvector u_vec = {n, x};
      dvector b_vec = {n, (REAL*)b};
      smoother_dcsr_gs(&u_vec, n-1, 0, -1, &mgl[lev].A, &b_vec, nu);
      break;
    }
    case SMOOTHER_GS_CF:
      cf_bwd_gs(&mgl[lev].A, aux->cf_order, n, b, x, nu);
      break;
    case SMOOTHER_JACOBI: {
      dvector u_vec = {n, x};
      dvector b_vec = {n, (REAL*)b};
      smoother_dcsr_jacobi(&u_vec, n-1, 0, -1, &mgl[lev].A, &b_vec, nu);
      break;
    }
    case SMOOTHER_L1DIAG:
    default:
      l1_smooth(&mgl[lev].A, aux->l1_diag, n, b, x, nu);
      break;
  }
}

/* ======================================================================
 *  Backslash recursive: pre-smooth + coarse correct
 * ====================================================================== */
static void backslash_rec(AMG_data* mgl, AMG_param* param, INT lev,
                          const REAL* b, REAL* x, INT nu) {
  INT nlev = mgl[0].num_levels;
  INT n = mgl[lev].A.row;

  if (lev == nlev - 1) {
    /* Coarsest level: ichol-preconditioned PCG */
    coarse_solve(&mgl[lev].A, b, x, &RS_AUX(mgl, lev)->L_ichol, n);
    return;
  }

  /* Pre-smoothing */
  smooth_fwd(mgl, param, lev, b, x, nu);

  /* Residual: r = b - A*x */
  REAL* r = (REAL*)calloc(n, sizeof(REAL));
  array_cp(n, (REAL*)b, r); dcsr_aAxpy(-1.0, &mgl[lev].A, x, r);

  /* Restrict: rc = P^T * r */
  INT nc = mgl[lev].P.col;
  REAL* rc = (REAL*)calloc(nc, sizeof(REAL));
  dcsr_mxv_trans(&mgl[lev].P, r, rc);

  /* Solve on coarse grid */
  REAL* ec = (REAL*)calloc(nc, sizeof(REAL));
  backslash_rec(mgl, param, lev + 1, rc, ec, nu);

  /* Prolongate: x = x + P * ec */
  dcsr_aAxpy(1.0, &mgl[lev].P, ec, x);

  free(r); free(rc); free(ec);
}

void rs_amg_backslash(AMG_data* mgl, AMG_param* param, INT lev,
                      const REAL* b, REAL* x, INT nu) {
  backslash_rec(mgl, param, lev, b, x, nu);
}

/* ======================================================================
 *  Fwdslash recursive: coarse correct + post-smooth
 * ====================================================================== */
static void fwdslash_rec(AMG_data* mgl, AMG_param* param, INT lev,
                         const REAL* b, REAL* x, INT nu) {
  INT nlev = mgl[0].num_levels;
  INT n = mgl[lev].A.row;

  if (lev == nlev - 1) {
    coarse_solve(&mgl[lev].A, b, x, &RS_AUX(mgl, lev)->L_ichol, n);
    return;
  }

  /* Residual: r = b - A*x */
  REAL* r = (REAL*)calloc(n, sizeof(REAL));
  array_cp(n, (REAL*)b, r); dcsr_aAxpy(-1.0, &mgl[lev].A, x, r);

  /* Restrict: rc = P^T * r */
  INT nc = mgl[lev].P.col;
  REAL* rc = (REAL*)calloc(nc, sizeof(REAL));
  dcsr_mxv_trans(&mgl[lev].P, r, rc);

  /* Solve on coarse grid */
  REAL* ec = (REAL*)calloc(nc, sizeof(REAL));
  fwdslash_rec(mgl, param, lev + 1, rc, ec, nu);

  /* Prolongate: x = x + P * ec */
  dcsr_aAxpy(1.0, &mgl[lev].P, ec, x);

  /* Post-smoothing */
  smooth_bwd(mgl, param, lev, b, x, nu);

  free(r); free(rc); free(ec);
}

void rs_amg_fwdslash(AMG_data* mgl, AMG_param* param, INT lev,
                     const REAL* b, REAL* x, INT nu) {
  fwdslash_rec(mgl, param, lev, b, x, nu);
}

/* ======================================================================
 *  Symmetric V-cycle: backslash + fwdslash
 * ====================================================================== */
void rs_amg_vcycle_precond(AMG_data* mgl, AMG_param* param,
                           const REAL* g, REAL* x) {
  INT n = mgl[0].A.row;
  INT nu = param->presmooth_iter;

  array_set(n, x, 0.0);
  rs_amg_backslash(mgl, param, 0, g, x, nu);

  /* r = g - A*x */
  REAL* r = (REAL*)calloc(n, sizeof(REAL));
  array_cp(n, (REAL*)g, r); dcsr_aAxpy(-1.0, &mgl[0].A, x, r);

  /* fwdslash on residual */
  REAL* corr = (REAL*)calloc(n, sizeof(REAL));
  rs_amg_fwdslash(mgl, param, 0, r, corr, nu);

  /* x = x + corr */
  array_axpy(n, 1.0, corr, x);

  free(r); free(corr);
}

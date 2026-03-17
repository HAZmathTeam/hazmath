/*
 * amg_wrapper.c - Shared library wrapper for AMG + ichol solver
 *
 * Authors: HAZmath (https://hazmath.net)
 * Created with the help of Claude (Anthropic)
 *
 * Exposes an opaque-handle API for FFI (Julia, Python, etc.).
 * Uses AMG_data and AMG_param from HAZmath.
 */
#include "amg_helping.h"

typedef struct {
  dCSRmat   A;
  AMG_data* mgl;
  AMG_param param;
  dCSRmat   L;
} amg_solver;

/* Preconditioner callback matching hazmath precond.fct signature */
typedef struct {
  AMG_data*      mgl;
  AMG_param*     param;
  const dCSRmat* A;
  const dCSRmat* L;
} precond_ctx;

static void precond_wrap(REAL* r, REAL* z, void* data) {
  const precond_ctx* ctx = (const precond_ctx*)data;
  rs_amg_ichol_precond(ctx->mgl, ctx->param, ctx->A, ctx->L, r, z);
}

/* Symmetrize a CSR matrix in-place: A = (A + A^T) / 2 */
static void symmetrize_dcsr(dCSRmat* A) {
  INT n = A->row;
  dCSRmat At;
  dcsr_alloc(n, n, A->nnz, &At);
  dcsr_transz(A, NULL, &At);

  dCSRmat Asum;
  dcsr_add(A, 0.5, &At, 0.5, &Asum);

  dcsr_free(A);
  dcsr_free(&At);
  *A = Asum;
}

/*
 * amg_solver_create - Build AMG hierarchy and ichol factor
 */
void* amg_solver_create(int n, int nnz,
                         const int* ia, const int* ja, const double* val,
                         double threshold, int min_size, int max_levels,
                         int nu) {
  amg_solver* s = (amg_solver*)calloc(1, sizeof(amg_solver));
  if (!s) return NULL;

  /* Copy CSR input into a dCSRmat */
  dcsr_alloc(n, n, nnz, &s->A);
  memcpy(s->A.IA, ia, (n + 1) * sizeof(INT));
  memcpy(s->A.JA, ja, nnz * sizeof(INT));
  memcpy(s->A.val, val, nnz * sizeof(REAL));

  /* Symmetrize: A = (A + A') / 2 */
  symmetrize_dcsr(&s->A);

  /* Set up AMG parameters */
  memset(&s->param, 0, sizeof(AMG_param));
  s->param.strong_coupled  = threshold;
  s->param.coarse_dof      = min_size;
  s->param.max_levels      = (max_levels > 0 && max_levels <= 20) ? max_levels : 20;
  s->param.smoother        = SMOOTHER_GS_CF;
  s->param.presmooth_iter  = (nu > 0) ? nu : 1;
  s->param.postsmooth_iter = (nu > 0) ? nu : 1;
  s->param.relaxation      = 0.5;

  /* Build AMG hierarchy */
  s->mgl = amg_data_create(s->param.max_levels);
  rs_amg_build_hierarchy(s->mgl, &s->param, &s->A);

  /* Incomplete Cholesky */
  ichol_compute(&s->A, &s->L);

  return (void*)s;
}

/*
 * amg_solver_update - Update solver for a spectrally equivalent matrix
 */
int amg_solver_update(void* handle, int n, int nnz,
                       const int* ia, const int* ja, const double* val) {
  amg_solver* s = (amg_solver*)handle;

  /* Build new CSR matrix */
  dCSRmat A_new;
  dcsr_alloc(n, n, nnz, &A_new);
  memcpy(A_new.IA, ia, (n + 1) * sizeof(INT));
  memcpy(A_new.JA, ja, nnz * sizeof(INT));
  memcpy(A_new.val, val, nnz * sizeof(REAL));

  /* Symmetrize */
  symmetrize_dcsr(&A_new);

  /* Free old fine-level ichol and matrix */
  dcsr_free(&s->L);
  dcsr_free(&s->A);
  s->A = A_new;

  /* Rebuild AMG hierarchy values */
  rs_amg_rebuild_values(s->mgl, &s->param, &s->A);

  /* Recompute fine-level ichol */
  ichol_compute(&s->A, &s->L);

  return 0;
}

/*
 * amg_solver_set_matrix - Replace solve matrix A and ichol, keep AMG hierarchy
 */
int amg_solver_set_matrix(void* handle, int n, int nnz,
                           const int* ia, const int* ja, const double* val) {
  amg_solver* s = (amg_solver*)handle;

  dCSRmat A_new;
  dcsr_alloc(n, n, nnz, &A_new);
  memcpy(A_new.IA, ia, (n + 1) * sizeof(INT));
  memcpy(A_new.JA, ja, nnz * sizeof(INT));
  memcpy(A_new.val, val, nnz * sizeof(REAL));

  symmetrize_dcsr(&A_new);

  dcsr_free(&s->L);
  dcsr_free(&s->A);
  s->A = A_new;

  ichol_compute(&s->A, &s->L);

  return 0;
}

/*
 * amg_solver_solve - Solve A*x = b using PCG with AMG+ichol
 */
int amg_solver_solve(void* handle,
                     const double* b, double* x,
                     double tol, int maxiter, int print_level) {
  amg_solver* s = (amg_solver*)handle;
  INT n = s->A.row;

  precond_ctx ctx;
  ctx.mgl   = s->mgl;
  ctx.param = &s->param;
  ctx.A     = &s->A;
  ctx.L     = &s->L;

  precond pc;
  pc.data = &ctx;
  pc.fct  = precond_wrap;

  dvector bv = {n, (REAL*)b};
  dvector uv = {n, x};

  return dcsr_pcg(&s->A, &bv, &uv, &pc, tol, maxiter, STOP_REL_PRECRES, print_level);
}

/*
 * amg_solver_set_smoother - Change smoother type and parameters
 */
void amg_solver_set_smoother(void* handle, int smoother, double omega) {
  amg_solver* s = (amg_solver*)handle;
  s->param.smoother = smoother;
  if (smoother == SMOOTHER_JACOBI)
    s->param.relaxation = omega;
}

/* Return matrix dimension */
int amg_solver_size(void* handle) {
  amg_solver* s = (amg_solver*)handle;
  return s->A.row;
}

/* Return number of AMG levels */
int amg_solver_num_levels(void* handle) {
  amg_solver* s = (amg_solver*)handle;
  return s->mgl[0].num_levels;
}

/* Free all resources */
void amg_solver_free(void* handle) {
  if (!handle) return;
  amg_solver* s = (amg_solver*)handle;
  dcsr_free(&s->L);
  rs_amg_free(s->mgl);
  dcsr_free(&s->A);
  free(s);
}

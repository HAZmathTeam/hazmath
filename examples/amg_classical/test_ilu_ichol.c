/*
 * test_ilu_ichol.c - Test incomplete factorization preconditioners
 *
 * Authors: HAZmath (https://hazmath.net)
 * Created with the help of Claude (Anthropic)
 *
 * Tests: ichol(0), icholt, ilu(0), ilut as PCG/PGMRES preconditioners.
 *
 * Usage: ./test_ilu_ichol <matrix_file> [tau]
 */
#include "hazmath.h"
#include <time.h>
#include <sys/resource.h>

/* --- ichol(0) preconditioner wrapper --- */
static void ichol0_precond_fct(REAL* r, REAL* z, void* data) {
  ichol_solve((const dCSRmat*)data, r, z);
}

/* --- icholt preconditioner wrapper --- */
static void icholt_precond_fct(REAL* r, REAL* z, void* data) {
  ichol_solve((const dCSRmat*)data, r, z);
}

/* --- ilu preconditioner data and wrapper --- */
typedef struct {
  const dCSRmat* L;
  const dCSRmat* U;
} ilu_precond_data;

static void ilu_precond_fct(REAL* r, REAL* z, void* data) {
  const ilu_precond_data* d = (const ilu_precond_data*)data;
  ilu_solve(d->L, d->U, r, z);
}

/* Load and symmetrize a matrix from a COO file */
static void load_and_symmetrize(const char* filename, dCSRmat* A) {
  dcoo_read_dcsr(filename, A);
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

int main(int argc, char** argv) {
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <matrix_file> [tau]\n", argv[0]);
    return 1;
  }

  const char* matfile = argv[1];
  REAL tau = (argc > 2) ? atof(argv[2]) : 0.01;

  /* Load matrix */
  fprintf(stderr, "Loading matrix from %s ...\n", matfile);
  clock_t t0 = clock();
  dCSRmat A;
  load_and_symmetrize(matfile, &A);
  INT n = A.row;
  double t_load = (double)(clock() - t0) / CLOCKS_PER_SEC;
  fprintf(stderr, "  n = %d, nnz = %d  (%.2f s)\n\n", n, A.nnz, t_load);

  /* Random RHS */
  REAL* b = (REAL*)malloc(n * sizeof(REAL));
  REAL* x = (REAL*)calloc(n, sizeof(REAL));
  dvector bv = {n, b};
  dvec_rand(n, &bv);
  dvector uv = {n, x};

  REAL* r = (REAL*)calloc(n, sizeof(REAL));
  REAL bnorm = sqrt(array_dotprod(n, b, b));
  REAL rnorm;
  REAL tol = 1e-6;
  INT maxiter = 500;

  /* ================================================================
   *  ichol(0)
   * ================================================================ */
  fprintf(stderr, "--- ichol(0) ---\n");
  dCSRmat L_ichol0;
  t0 = clock();
  ichol_compute(&A, &L_ichol0);
  double t_setup = (double)(clock() - t0) / CLOCKS_PER_SEC;
  fprintf(stderr, "  setup: %.2f s, nnz(L) = %d\n", t_setup, L_ichol0.nnz);

  precond pc;
  pc.data = &L_ichol0;
  pc.fct  = ichol0_precond_fct;

  array_set(n, x, 0.0);
  t0 = clock();
  INT iter = dcsr_pcg(&A, &bv, &uv, &pc, tol, maxiter, STOP_REL_PRECRES, 1);
  double t_solve = (double)(clock() - t0) / CLOCKS_PER_SEC;
  array_cp(n, b, r); dcsr_aAxpy(-1.0, &A, x, r);
  rnorm = sqrt(array_dotprod(n, r, r));
  fprintf(stderr, "  PCG: iter = %d, relres = %.2e, solve = %.2f s\n\n",
          iter, rnorm / bnorm, t_solve);

  /* ================================================================
   *  icholt(tau)
   * ================================================================ */
  fprintf(stderr, "--- icholt(tau=%.4f) ---\n", tau);
  dCSRmat L_icholt;
  t0 = clock();
  icholt_compute(&A, tau, &L_icholt);
  t_setup = (double)(clock() - t0) / CLOCKS_PER_SEC;
  fprintf(stderr, "  setup: %.2f s, nnz(L) = %d\n", t_setup, L_icholt.nnz);

  precond pc_icholt;
  pc_icholt.data = &L_icholt;
  pc_icholt.fct  = icholt_precond_fct;

  array_set(n, x, 0.0);
  t0 = clock();
  iter = dcsr_pcg(&A, &bv, &uv, &pc_icholt, tol, maxiter, STOP_REL_PRECRES, 1);
  t_solve = (double)(clock() - t0) / CLOCKS_PER_SEC;
  array_cp(n, b, r); dcsr_aAxpy(-1.0, &A, x, r);
  rnorm = sqrt(array_dotprod(n, r, r));
  fprintf(stderr, "  PCG: iter = %d, relres = %.2e, solve = %.2f s\n\n",
          iter, rnorm / bnorm, t_solve);

  /* ================================================================
   *  ILU(0)
   * ================================================================ */
  fprintf(stderr, "--- ILU(0) ---\n");
  dCSRmat L_ilu0, U_ilu0;
  t0 = clock();
  ilu_compute(&A, &L_ilu0, &U_ilu0);
  t_setup = (double)(clock() - t0) / CLOCKS_PER_SEC;
  fprintf(stderr, "  setup: %.2f s, nnz(L) = %d, nnz(U) = %d\n",
          t_setup, L_ilu0.nnz, U_ilu0.nnz);

  ilu_precond_data ilu0_data = {&L_ilu0, &U_ilu0};
  precond pc_ilu0;
  pc_ilu0.data = &ilu0_data;
  pc_ilu0.fct  = ilu_precond_fct;

  array_set(n, x, 0.0);
  t0 = clock();
  iter = dcsr_pcg(&A, &bv, &uv, &pc_ilu0, tol, maxiter, STOP_REL_PRECRES, 1);
  t_solve = (double)(clock() - t0) / CLOCKS_PER_SEC;
  array_cp(n, b, r); dcsr_aAxpy(-1.0, &A, x, r);
  rnorm = sqrt(array_dotprod(n, r, r));
  fprintf(stderr, "  PCG: iter = %d, relres = %.2e, solve = %.2f s\n\n",
          iter, rnorm / bnorm, t_solve);

  /* ================================================================
   *  ILUT(tau)
   * ================================================================ */
  fprintf(stderr, "--- ILUT(tau=%.4f) ---\n", tau);
  dCSRmat L_ilut, U_ilut;
  t0 = clock();
  ilut_compute(&A, tau, &L_ilut, &U_ilut);
  t_setup = (double)(clock() - t0) / CLOCKS_PER_SEC;
  fprintf(stderr, "  setup: %.2f s, nnz(L) = %d, nnz(U) = %d\n",
          t_setup, L_ilut.nnz, U_ilut.nnz);

  ilu_precond_data ilut_data = {&L_ilut, &U_ilut};
  precond pc_ilut;
  pc_ilut.data = &ilut_data;
  pc_ilut.fct  = ilu_precond_fct;

  array_set(n, x, 0.0);
  t0 = clock();
  iter = dcsr_pcg(&A, &bv, &uv, &pc_ilut, tol, maxiter, STOP_REL_PRECRES, 1);
  t_solve = (double)(clock() - t0) / CLOCKS_PER_SEC;
  array_cp(n, b, r); dcsr_aAxpy(-1.0, &A, x, r);
  rnorm = sqrt(array_dotprod(n, r, r));
  fprintf(stderr, "  PCG: iter = %d, relres = %.2e, solve = %.2f s\n\n",
          iter, rnorm / bnorm, t_solve);

  /* ================================================================
   *  Summary
   * ================================================================ */
  fprintf(stderr, "========================================\n");
  fprintf(stderr, "  n = %d, tau = %.4f\n", n, tau);
  fprintf(stderr, "========================================\n");

  /* Cleanup */
  dcsr_free(&L_ichol0);
  dcsr_free(&L_icholt);
  dcsr_free(&L_ilu0); dcsr_free(&U_ilu0);
  dcsr_free(&L_ilut); dcsr_free(&U_ilut);
  dcsr_free(&A);
  free(b); free(x); free(r);

  /* Peak RAM usage */
  {
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
#ifdef __APPLE__
    double peak_mb = (double)ru.ru_maxrss / (1024.0 * 1024.0); /* bytes on macOS */
#else
    double peak_mb = (double)ru.ru_maxrss / 1024.0; /* KB on Linux */
#endif
    fprintf(stderr, "\nPeak RSS = %.2f MB\n", peak_mb);
  }
  return 0;
}

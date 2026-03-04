/*
 * test_amg.c - Test AMG V-cycle and AMG+ichol preconditioners
 *
 * Authors: HAZmath (https://hazmath.net)
 * Created with the help of Claude (Anthropic)
 *
 * Usage: ./test_amg <solve_matrix> [amg_matrix] [threshold] [nu] [min_size]
 *
 * Uses AMG_data and AMG_param from HAZmath.
 */
#include "amg_helping.h"

/* Preconditioner wrappers for hazmath precond struct */

typedef struct {
  AMG_data *mgl;
  AMG_param *param;
  const dCSRmat* A;
  const dCSRmat* L;
} pcg_precond_data;

static void vcycle_precond_wrap(REAL* r, REAL* z, void* data) {
  const pcg_precond_data* d = (const pcg_precond_data*)data;
  rs_amg_vcycle_precond(d->mgl, d->param, r, z);
}

static void ichol_precond_wrap(REAL* r, REAL* z, void* data) {
  const pcg_precond_data* d = (const pcg_precond_data*)data;
  rs_amg_ichol_precond(d->mgl, d->param, d->A, d->L, r, z);
}

/* Load and symmetrize a matrix from a COO file */
static int load_matrix(const char* filename, dCSRmat* A) {
  dcoo_read_dcsr(filename, A);

  /* Symmetrize: A = (A + A') / 2 */
  INT n = A->row;
  dCSRmat At;
  dcsr_alloc(n, n, A->nnz, &At);
  dcsr_transz(A, NULL, &At);

  dCSRmat Asum;
  dcsr_add(A, 0.5, &At, 0.5, &Asum);

  dcsr_free(A);
  dcsr_free(&At);
  *A = Asum;

  return 0;
}

int main(int argc, char** argv) {
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <solve_matrix> [amg_matrix] [threshold] [nu] [min_size]\n", argv[0]);
    return 1;
  }

  /* Parse arguments: detect whether arg2 is a filename or a number */
  const char* solve_file = argv[1];
  const char* amg_file   = NULL;
  int arg_offset = 2;

  if (argc >= 3) {
    char c = argv[2][0];
    if (c != '.' && (c < '0' || c > '9')) {
      amg_file = argv[2];
      arg_offset = 3;
    }
  }

  REAL threshold = (argc > arg_offset)     ? atof(argv[arg_offset])     : 0.25;
  INT  nu        = (argc > arg_offset + 1) ? atoi(argv[arg_offset + 1]) : 2;
  INT  min_size  = (argc > arg_offset + 2) ? atoi(argv[arg_offset + 2]) : 4096;

  int two_matrix = (amg_file != NULL);

  /* Load solve matrix A */
  fprintf(stderr, "Loading solve matrix from %s ...\n", solve_file);
  clock_t t0 = clock();

  dCSRmat A;
  if (load_matrix(solve_file, &A) != 0) return 1;

  INT n = A.row;
  double t_load = (double)(clock() - t0) / CLOCKS_PER_SEC;
  fprintf(stderr, "  n = %d, nnz = %d  (%.2f s)\n", n, A.nnz, t_load);

  /* Load AMG matrix Ap (or alias to A) */
  dCSRmat Ap_storage;
  dCSRmat* Ap;
  if (two_matrix) {
    fprintf(stderr, "Loading AMG matrix from %s ...\n", amg_file);
    t0 = clock();
    if (load_matrix(amg_file, &Ap_storage) != 0) return 1;
    double t_load2 = (double)(clock() - t0) / CLOCKS_PER_SEC;
    fprintf(stderr, "  n = %d, nnz = %d  (%.2f s)\n",
            Ap_storage.row, Ap_storage.nnz, t_load2);
    if (Ap_storage.row != n) {
      fprintf(stderr, "ERROR: matrix dimensions do not match (%d vs %d)\n",
              Ap_storage.row, n);
      return 1;
    }
    Ap = &Ap_storage;
  } else {
    Ap = &A;
  }

  /* RHS: try rhs.txt, fall back to random */
  REAL* b = (REAL*)malloc(n * sizeof(REAL));
  REAL* x = (REAL*)calloc(n, sizeof(REAL));
  {
    FILE* frhs = fopen("rhs.txt", "r");
    if (frhs) {
      INT nrhs;
      if (fscanf(frhs, "%d", &nrhs) == 1 && nrhs == n) {
        for (INT i = 0; i < n; i++) {
          if (fscanf(frhs, "%lf", &b[i]) != 1) {
            fprintf(stderr, "Error reading rhs.txt at entry %d, using random RHS\n", i);
            dvector bv = {n, b};
            dvec_rand(n, &bv);
            break;
          }
        }
        fprintf(stderr, "RHS loaded from rhs.txt\n");
      } else {
        fprintf(stderr, "rhs.txt dimension mismatch (expected %d), using random RHS\n", n);
        dvector bv = {n, b};
        dvec_rand(n, &bv);
      }
      fclose(frhs);
    } else {
      fprintf(stderr, "rhs.txt not found, using random RHS\n");
      dvector bv = {n, b};
      dvec_rand(n, &bv);
    }
  }

  /* Set up AMG parameters */
  AMG_param param;
  memset(&param, 0, sizeof(param));
  param.strong_coupled  = threshold;
  param.coarse_dof      = min_size;
  param.max_levels      = 20;
  param.smoother        = SMOOTHER_GS_CF;
  param.presmooth_iter  = nu;
  param.postsmooth_iter = nu;
  param.relaxation      = 0.5;

  /* Build AMG hierarchy from Ap */
  AMG_data *mgl = amg_data_create(param.max_levels);

  if (two_matrix)
    fprintf(stderr, "\nBuilding AMG hierarchy from %s ...\n", amg_file);

  t0 = clock();
  rs_amg_build_hierarchy(mgl, &param, Ap);
  double t_hier = (double)(clock() - t0) / CLOCKS_PER_SEC;
  fprintf(stderr, "Hierarchy built in %.2f s\n", t_hier);

  /* Report nnz per F-row in P at each level */
  for (INT lev = 0; lev < mgl[0].num_levels - 1; lev++) {
    const dCSRmat* P = &mgl[lev].P;
    const INT* mis = RS_AUX(mgl, lev)->mis;
    INT nf = 0, min_nnz = P->row, max_nnz = 0;
    long total_nnz = 0;
    for (INT i = 0; i < P->row; i++) {
      if (mis[i]) continue;
      nf++;
      INT rnnz = P->IA[i + 1] - P->IA[i];
      total_nnz += rnnz;
      if (rnnz < min_nnz) min_nnz = rnnz;
      if (rnnz > max_nnz) max_nnz = rnnz;
    }
    fprintf(stderr, "  Level %d P: %d F-rows, nnz/F-row: min=%d, max=%d, avg=%.1f\n",
            lev + 1, nf, min_nnz, max_nnz, nf > 0 ? (double)total_nnz / nf : 0.0);
  }

  /* Incomplete Cholesky */
  dCSRmat L;
  t0 = clock();
  if (two_matrix) {
    fprintf(stderr, "Computing ichol of AMG matrix %s ...\n", amg_file);
    ichol_compute(Ap, &L);
  } else {
    ichol_compute(&A, &L);
  }
  double t_ichol = (double)(clock() - t0) / CLOCKS_PER_SEC;
  fprintf(stderr, "ichol computed in %.2f s (nnz(L) = %d)\n", t_ichol, L.nnz);

  REAL tol = 1e-6;
  INT maxiter = 200;

  /* Preconditioner data */
  pcg_precond_data pdata;
  pdata.mgl   = mgl;
  pdata.param = &param;
  pdata.A     = two_matrix ? Ap : &A;
  pdata.L     = &L;

  /* dvector wrappers for hazmath dcsr_pcg */
  dvector bv = {n, b};
  dvector uv = {n, x};
  REAL* r = (REAL*)calloc(n, sizeof(REAL));
  REAL bnorm = sqrt(array_dotprod(n, b, b));
  REAL rnorm;

  /* --- PCG with AMG+ichol --- */
  if (two_matrix)
    fprintf(stderr, "\nPCG solve of %s with AMG(%s)+ichol(%s):\n",
            solve_file, amg_file, amg_file);
  else
    fprintf(stderr, "\nPCG with AMG+ichol preconditioner:\n");

  array_set(n, x, 0.0);
  precond pc_ichol;
  pc_ichol.data = &pdata;
  pc_ichol.fct  = ichol_precond_wrap;
  t0 = clock();
  INT iter1 = dcsr_pcg(&A, &bv, &uv, &pc_ichol, tol, maxiter, STOP_REL_PRECRES, 1);
  double t1 = (double)(clock() - t0) / CLOCKS_PER_SEC;

  array_cp(n, b, r); dcsr_aAxpy(-1.0, &A, x, r);
  rnorm = sqrt(array_dotprod(n, r, r));
  fprintf(stderr, "  iterations = %d, relres = %.2e, time = %.2f s\n",
          iter1, rnorm / bnorm, t1);

  /* --- PCG with AMG V-cycle only --- */
  if (two_matrix)
    fprintf(stderr, "\nPCG solve of %s with AMG(%s) V-cycle only:\n",
            solve_file, amg_file);
  else
    fprintf(stderr, "\nPCG with AMG V-cycle preconditioner:\n");

  array_set(n, x, 0.0);
  precond pc_vcycle;
  pc_vcycle.data = &pdata;
  pc_vcycle.fct  = vcycle_precond_wrap;
  t0 = clock();
  INT iter2 = dcsr_pcg(&A, &bv, &uv, &pc_vcycle, tol, maxiter, STOP_REL_PRECRES, 1);
  double t2 = (double)(clock() - t0) / CLOCKS_PER_SEC;

  array_cp(n, b, r); dcsr_aAxpy(-1.0, &A, x, r);
  rnorm = sqrt(array_dotprod(n, r, r));
  fprintf(stderr, "  iterations = %d, relres = %.2e, time = %.2f s\n",
          iter2, rnorm / bnorm, t2);

  /* Summary */
  fprintf(stderr, "\n========================================\n");
  if (two_matrix)
    fprintf(stderr, "  AMG from %s, solve %s\n", amg_file, solve_file);
  fprintf(stderr, "  %-25s %6s %12s %12s\n", "Preconditioner", "Iter", "Setup(s)", "Solve(s)");
  fprintf(stderr, "  %-25s %6d %12.2f %12.2f\n", "AMG+ichol", iter1, t_hier + t_ichol, t1);
  fprintf(stderr, "  %-25s %6d %12.2f %12.2f\n", "AMG V-cycle only", iter2, t_hier, t2);
  fprintf(stderr, "========================================\n");

  /* Cleanup */
  if (two_matrix) dcsr_free(&Ap_storage);
  free(b); free(x); free(r);
  dcsr_free(&L);
  rs_amg_free(mgl);
  dcsr_free(&A);

  return 0;
}

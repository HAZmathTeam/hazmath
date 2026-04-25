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
#include <string.h>
#include <sys/resource.h>

/* Preconditioner wrappers for hazmath precond struct */

typedef struct {
  AMG_data* mgl;
  AMG_param* param;
  const dCSRmat* A;
  const dCSRmat* L;
} pcg_precond_data;

/* RS AMG V-cycle preconditioner */
static void rs_vcycle_precond_wrap(REAL* r, REAL* z, void* data) {
  const pcg_precond_data* d = (const pcg_precond_data*)data;
  rs_amg_vcycle_precond(d->mgl, d->param, r, z);
}

/* RS AMG + ichol preconditioner */
static void ichol_precond_wrap(REAL* r, REAL* z, void* data) {
  const pcg_precond_data* d = (const pcg_precond_data*)data;
  rs_amg_ichol_precond(d->mgl, d->param, d->A, d->L, r, z);
}

/* UA AMG W-cycle preconditioner (uses hazmath mgcycle) */
static void ua_wcycle_precond_wrap(REAL* r, REAL* z, void* data) {
  pcg_precond_data* d = (pcg_precond_data*)data;
  AMG_data* mgl = d->mgl;
  INT m = mgl[0].A.row;
  mgl->b.row = m; array_cp(m, r, mgl->b.val);
  mgl->x.row = m; dvec_set(m, &mgl->x, 0.0);
  mgcycle(mgl, d->param);
  array_cp(m, mgl->x.val, z);
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
    fprintf(stderr, "Usage: %s [--cycle=V|W|UA] <solve_matrix> [amg_matrix] [threshold] [nu] [min_size]\n", argv[0]);
    fprintf(stderr, "  --cycle=V   RS AMG, V-cycle (default)\n");
    fprintf(stderr, "  --cycle=W   RS AMG, W-cycle\n");
    fprintf(stderr, "  --cycle=UA  UA AMG, W-cycle\n");
    return 1;
  }

  /* Parse --cycle option */
  int cycle_mode = 0;  /* 0=RS V-cycle, 1=RS W-cycle, 2=UA W-cycle */
  int iarg = 1;
  if (argc > 1 && strncmp(argv[1], "--cycle=", 8) == 0) {
    const char* val = argv[1] + 8;
    if (strcmp(val, "W") == 0 || strcmp(val, "w") == 0)
      cycle_mode = 1;
    else if (strcmp(val, "UA") == 0 || strcmp(val, "ua") == 0)
      cycle_mode = 2;
    else if (strcmp(val, "V") == 0 || strcmp(val, "v") == 0)
      cycle_mode = 0;
    else {
      fprintf(stderr, "Unknown cycle type '%s'. Use V, W, or UA.\n", val);
      return 1;
    }
    iarg = 2;
  }

  if (iarg >= argc) {
    fprintf(stderr, "Missing matrix file. Run with no arguments for usage.\n");
    return 1;
  }

  /* Parse remaining arguments */
  const char* solve_file = argv[iarg];
  const char* amg_file   = NULL;
  int arg_offset = iarg + 1;

  if (arg_offset < argc) {
    char c = argv[arg_offset][0];
    if (c != '.' && (c < '0' || c > '9') && c != '-') {
      amg_file = argv[arg_offset];
      arg_offset++;
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

  int use_ua = (cycle_mode == 2);

  AMG_param param;
  AMG_data* mgl;
  double t_hier;
  const char* amg_label;

  if (use_ua) {
    /* ---- UA AMG with W-cycle ---- */
    param_amg_init(&param);
    param.cycle_type      = W_CYCLE;
    param.max_levels      = 20;
    param.coarse_dof      = (min_size < 200) ? min_size : 200;
    param.print_level     = PRINT_MORE;

    mgl = amg_data_create(param.max_levels);
    dcsr_alloc(Ap->row, Ap->col, Ap->nnz, &mgl[0].A);
    dcsr_cp((dCSRmat*)Ap, &mgl[0].A);
    mgl[0].b = dvec_create(n);
    mgl[0].x = dvec_create(n);

    fprintf(stderr, "\nBuilding UA AMG hierarchy ...\n");
    t0 = clock();
    amg_setup_ua(mgl, &param);
    t_hier = (double)(clock() - t0) / CLOCKS_PER_SEC;
    fprintf(stderr, "UA hierarchy built in %.2f s (%d levels)\n",
            t_hier, mgl[0].num_levels);
    amg_label = "UA W-cycle";
  } else {
    /* ---- RS AMG with V or W-cycle ---- */
    memset(&param, 0, sizeof(param));
    param.strong_coupled  = threshold;
    param.coarse_dof      = min_size;
    param.max_levels      = 20;
    param.smoother        = SMOOTHER_GS_CF;
    param.presmooth_iter  = nu;
    param.postsmooth_iter = nu;
    param.relaxation      = 0.5;
    param.cycle_type      = (cycle_mode == 1) ? 2 : 1;  /* 1=V-cycle, 2=W-cycle */

    mgl = amg_data_create(param.max_levels);

    if (two_matrix)
      fprintf(stderr, "\nBuilding RS AMG hierarchy from %s ...\n", amg_file);

    t0 = clock();
    rs_amg_build_hierarchy(mgl, &param, Ap);
    t_hier = (double)(clock() - t0) / CLOCKS_PER_SEC;
    fprintf(stderr, "RS hierarchy built in %.2f s (%d levels)\n",
            t_hier, mgl[0].num_levels);

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
    amg_label = (param.cycle_type >= 2) ? "RS W-cycle" : "RS V-cycle";
  }

  REAL tol = 1e-6;
  INT maxiter = 200;

  /* Preconditioner data */
  pcg_precond_data pdata;
  pdata.mgl   = mgl;
  pdata.param = &param;
  pdata.A     = two_matrix ? Ap : &A;
  pdata.L     = NULL;

  /* dvector wrappers for hazmath dcsr_pcg */
  dvector bv = {n, b};
  dvector uv = {n, x};
  REAL* r = (REAL*)calloc(n, sizeof(REAL));
  REAL bnorm = sqrt(array_dotprod(n, b, b));
  REAL rnorm;

  /* --- PCG with AMG preconditioner --- */
  fprintf(stderr, "\nPCG with %s preconditioner:\n", amg_label);

  array_set(n, x, 0.0);
  precond pc;
  pc.data = &pdata;
  pc.fct  = use_ua ? ua_wcycle_precond_wrap : rs_vcycle_precond_wrap;
  t0 = clock();
  INT iter = dcsr_pcg(&A, &bv, &uv, &pc, tol, maxiter, STOP_REL_PRECRES, 1);
  double t_solve = (double)(clock() - t0) / CLOCKS_PER_SEC;

  array_cp(n, b, r); dcsr_aAxpy(-1.0, &A, x, r);
  rnorm = sqrt(array_dotprod(n, r, r));
  fprintf(stderr, "  iterations = %d, relres = %.2e, time = %.2f s\n",
          iter, rnorm / bnorm, t_solve);

  /* Summary */
  fprintf(stderr, "\n========================================\n");
  if (two_matrix)
    fprintf(stderr, "  AMG from %s, solve %s\n", amg_file, solve_file);
  fprintf(stderr, "  %-25s %6s %12s %12s\n", "Preconditioner", "Iter", "Setup(s)", "Solve(s)");
  fprintf(stderr, "  %-25s %6d %12.2f %12.2f\n", amg_label, iter, t_hier, t_solve);
  fprintf(stderr, "========================================\n");

  /* Cleanup */
  if (two_matrix) dcsr_free(&Ap_storage);
  free(b); free(x); free(r);
  if (use_ua)
    amg_data_free(mgl, &param);
  else
    rs_amg_free(mgl);
  dcsr_free(&A);

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

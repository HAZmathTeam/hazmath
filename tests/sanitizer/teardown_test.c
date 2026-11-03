/* Build each HX variant and an RA preconditioner, apply it (the 2D
 * multiplicative HX-div five times), and free everything with
 * free_precond_*. Under ASan/UBSan/LSan nothing may be reported. Also: an
 * HX setup with the AMLI cycle returns NULL, and param_amg_cp copies the
 * AMLI coefficients. Run by tests/run_all_tests.sh. */
#include "poisson.h"
precond* create_precond_hxcurl(dCSRmat*, dCSRmat*, dCSRmat*, SHORT, AMG_param*);
precond* create_precond_hxdiv_3D(dCSRmat*, dCSRmat*, dCSRmat*, dCSRmat*, SHORT,
                                 AMG_param*);
precond* create_precond_hxdiv_2D(dCSRmat*, dCSRmat*, dCSRmat*, SHORT,
                                 AMG_param*);
precond* create_precond_ra(dCSRmat*, dCSRmat*, REAL, REAL, REAL, REAL, REAL,
                           REAL, REAL, AMG_param*);
void free_precond_hxcurl(precond*);
void free_precond_hxdiv(precond*);
void free_precond_ra(precond*);

static dCSRmat identity(INT n)
{
  dCSRmat I = dcsr_create(n, n, n);
  for (INT i = 0; i < n; ++i)
  {
    I.IA[i] = i;
    I.JA[i] = i;
    I.val[i] = 1.0;
  }
  I.IA[n] = n;
  return I;
}

static dCSRmat tridiag(INT n, REAL diag, REAL off)
{
  dCSRmat A = dcsr_create(n, n, 3 * n - 2);
  INT k = 0;
  A.IA[0] = 0;
  for (INT i = 0; i < n; ++i)
  {
    if (i > 0) { A.JA[k] = i - 1; A.val[k++] = off; }
    A.JA[k] = i; A.val[k++] = diag;
    if (i < n - 1) { A.JA[k] = i + 1; A.val[k++] = off; }
    A.IA[i + 1] = k;
  }
  return A;
}

static REAL apply(precond* pc, INT n)
{
  REAL* r = calloc(n, sizeof(REAL));
  REAL* z = calloc(n, sizeof(REAL));
  for (INT i = 0; i < n; ++i) r[i] = 1.0 + (i % 7);
  pc->fct(r, z, pc->data);
  REAL s = 0;
  for (INT i = 0; i < n; ++i) s += z[i] * z[i];
  free(r);
  free(z);
  return sqrt(s);
}

int main(void)
{
  dCSRmat A = poisson_2d(16);
  const INT n = A.row;
  dCSRmat I1 = identity(n), I2 = identity(n), I3 = identity(n);
  const SHORT curl[2] = {PREC_HX_CURL_A, PREC_HX_CURL_M};
  const SHORT div[2] = {PREC_HX_DIV_A, PREC_HX_DIV_M};
  AMG_param param;

  for (int k = 0; k < 2; ++k)
  {
    param_amg_init(&param);
    precond* pc = create_precond_hxcurl(&A, &I1, &I2, curl[k], &param);
    printf("HX-curl %s: |B r| = %.6e\n", k ? "mult" : "add ", apply(pc, n));
    free_precond_hxcurl(pc);

    param_amg_init(&param);
    pc = create_precond_hxdiv_3D(&A, &I1, &I2, &I3, div[k], &param);
    printf("HX-div 3D %s: |B r| = %.6e\n", k ? "mult" : "add ", apply(pc, n));
    free_precond_hxdiv(pc);
  }
  /* 2D, both variants; the multiplicative apply leaked on every call
   * (issue 10), so apply it several times */
  precond* pc = NULL;
  for (int k = 0; k < 2; ++k)
  {
    param_amg_init(&param);
    pc = create_precond_hxdiv_2D(&A, &I1, &I2, div[k], &param);
    REAL nrm = 0;
    for (int it = 0; it < 5; ++it) nrm = apply(pc, n);
    printf("HX-div 2D %s: |B r| = %.6e\n", k ? "mult" : "add ", nrm);
    free_precond_hxdiv(pc);
  }

  /* AMLI is rejected by the HX factories */
  param_amg_init(&param);
  param.cycle_type = AMLI_CYCLE;
  pc = create_precond_hxcurl(&A, &I1, &I2, PREC_HX_CURL_A, &param);
  printf("HX-curl with AMLI returns %s\n", pc ? "a preconditioner" : "NULL");

  /* param_amg_cp allocates the destination's AMLI coefficients */
  AMG_param p1, p2;
  param_amg_init(&p1);
  param_amg_init(&p2);
  p1.amli_degree = 2;
  p1.amli_coef = malloc(3 * sizeof(REAL));
  p1.amli_coef[0] = 1; p1.amli_coef[1] = 2; p1.amli_coef[2] = 3;
  param_amg_cp(&p1, &p2);
  printf("param_amg_cp: %s copy, coef = %g %g %g\n",
         p2.amli_coef != p1.amli_coef ? "separate" : "SHARED",
         p2.amli_coef[0], p2.amli_coef[1], p2.amli_coef[2]);
  free(p1.amli_coef);
  free(p2.amli_coef);

  /* RA, s = -0.5; its setup still leaks (issue 18, step 2) */
  const INT m = 63;
  const REAL h = 1.0 / (m + 1);
  dCSRmat K = tridiag(m, 2.0 / h, -1.0 / h);
  dCSRmat M = tridiag(m, 4.0 * h / 6, h / 6);
  param_amg_init(&param);
  pc = create_precond_ra(&K, &M, -0.5, 0.0, 1.0, 0.0, h * h / (1.1 * 12.0),
                         1.0, 1e-10, &param);
  printf("RA: |B r| = %.6e\n", apply(pc, m));
  free_precond_ra(pc);

  dcsr_free(&A);
  dcsr_free(&I1);
  dcsr_free(&I2);
  dcsr_free(&I3);
  dcsr_free(&K);
  dcsr_free(&M);
  return 0;
}

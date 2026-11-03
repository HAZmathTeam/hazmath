/* Failed HX and RA setups free what they built: max_levels = 0 makes the
 * first AMG setup in each factory fail (amg_setup_ua returns
 * ERROR_INPUT_PAR). Under LSan no leaks should be reported. */
#include "poisson.h"
precond* create_precond_hxcurl(dCSRmat*, dCSRmat*, dCSRmat*, SHORT, AMG_param*);
precond* create_precond_hxdiv_3D(dCSRmat*, dCSRmat*, dCSRmat*, dCSRmat*, SHORT,
                                 AMG_param*);
precond* create_precond_hxdiv_2D(dCSRmat*, dCSRmat*, dCSRmat*, SHORT,
                                 AMG_param*);
precond* create_precond_ra(dCSRmat*, dCSRmat*, REAL, REAL, REAL, REAL, REAL,
                           REAL, REAL, AMG_param*);

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

int main(void)
{
  dCSRmat A = poisson_2d(8);
  const INT n = A.row;
  dCSRmat I1 = identity(n), I2 = identity(n), I3 = identity(n);
  AMG_param param;

  param_amg_init(&param);
  param.max_levels = 0;
  printf("HX-curl:   %s\n", create_precond_hxcurl(&A, &I1, &I2, PREC_HX_CURL_A, &param) ? "built" : "NULL");
  param_amg_init(&param);
  param.max_levels = 0;
  printf("HX-div 3D: %s\n", create_precond_hxdiv_3D(&A, &I1, &I2, &I3, PREC_HX_DIV_A, &param) ? "built" : "NULL");
  param_amg_init(&param);
  param.max_levels = 0;
  printf("HX-div 2D: %s\n", create_precond_hxdiv_2D(&A, &I1, &I2, PREC_HX_DIV_A, &param) ? "built" : "NULL");

  const INT m = 63;
  const REAL h = 1.0 / (m + 1);
  dCSRmat K = tridiag(m, 2.0 / h, -1.0 / h);
  dCSRmat M = tridiag(m, 4.0 * h / 6, h / 6);
  param_amg_init(&param);
  param.max_levels = 0;
  printf("RA:        %s\n", create_precond_ra(&K, &M, -0.5, 0.0, 1.0, 0.0, h * h / (1.1 * 12.0),
                                              1.0, 1e-10, &param) ? "built" : "NULL");

  dcsr_free(&A);
  dcsr_free(&I1);
  dcsr_free(&I2);
  dcsr_free(&I3);
  dcsr_free(&K);
  dcsr_free(&M);
  return 0;
}

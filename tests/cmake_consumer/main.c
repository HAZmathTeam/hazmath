/*
 * Solve a 1D Laplacian with HAZmath AMG and check the relative residual.
 * Used by CI to check that an installed HAZmath is usable.
 */
#include "hazmath.h"

int main(void)
{
  const INT n = 1000;
  dCSRmat A = dcsr_create(n, n, 3 * n - 2);
  INT k = 0;
  for (INT i = 0; i < n; ++i)
  {
    A.IA[i] = k;
    if (i > 0)
    {
      A.JA[k] = i - 1;
      A.val[k++] = -1.0;
    }
    A.JA[k] = i;
    A.val[k++] = 2.0;
    if (i < n - 1)
    {
      A.JA[k] = i + 1;
      A.val[k++] = -1.0;
    }
  }
  A.IA[n] = k;

  dvector b = dvec_create(n), x = dvec_create(n), r = dvec_create(n);
  dvec_set(n, &b, 1.0);
  dvec_set(n, &x, 0.0);

  AMG_param amgparam;
  param_amg_init(&amgparam);
  amgparam.print_level = 0;
  amgparam.maxit = 500;
  amgparam.tol = 1e-10;
  linear_solver_amg(&A, &b, &x, &amgparam);

  /* r = b - A x */
  dvec_cp(&b, &r);
  dcsr_aAxpy(-1.0, &A, x.val, r.val);
  const REAL relres = dvec_norm2(&r) / dvec_norm2(&b);
  printf("relative residual = %e\n", relres);

  dcsr_free(&A);
  dvec_free(&b);
  dvec_free(&x);
  dvec_free(&r);
  return relres < 1e-8 ? 0 : 1;
}

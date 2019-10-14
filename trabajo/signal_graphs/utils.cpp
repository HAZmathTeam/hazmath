#include "hazmath_include.h"

REAL *initializeRhs(dCSRmat *A, int num_iterations) {
  assert(A->row == A->col);
  int n = A->row;
  dvector *f = (dvector *)malloc(sizeof(dvector)),
          *zero = (dvector *)malloc(sizeof(dvector));
  dvec_alloc(n, f);
  dvec_alloc(n, zero);
  dvec_set(n, zero, 0.0);

  for (int i = 0; i < n; ++i) {
    f->val[i] = (double)i / (n - 1);
  }

  dvec_rand(n, f);
  REAL sum = 0.0;
  for (int i = 0; i < n; ++i) {
    sum += f->val[i];
  }
  for (int i = 0; i < n; ++i) {
    f->val[i] -= sum / n;
  }
  dvec_ax(1.0 / dvec_norm2(f), f);

  for (int i = 0; i < num_iterations; ++i) {
    smoother_dcsr_sgs(f, A, zero, 1);
    REAL sum = 0;
    for (int i = 0; i < n; ++i) {
      sum += f->val[i];
    }
    for (int i = 0; i < n; ++i) {
      f->val[i] -= sum / n;
    }
    dvec_ax(1.0 / dvec_norm2(f), f);
    // dvector *g = (dvector *)malloc(sizeof(dvector));
    // dvec_alloc(n, g);
    // dcsr_mxv(A, f->val, g->val);
    // cout << "lambda: " << dvec_dotprod(f, g) / dvec_dotprod(f, f) << endl;
  }
  free(zero);
  return f->val;
}

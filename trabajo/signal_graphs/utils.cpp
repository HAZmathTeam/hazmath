#include "hazmath_include.h"
#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

REAL *initializeRhs(dCSRmat *L, int num_iterations) {
  assert(L->row == L->col);
  int n = L->row;
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
    smoother_dcsr_sgs(f, L, zero, 1);
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
    // dcsr_mxv(L, f->val, g->val);
    // cout << "lambda: " << dvec_dotprod(f, g) / dvec_dotprod(f, f) << endl;
  }
  free(zero);
  return f->val;
}

extern "C" {
void dsyevr_(char *JOBZ, char *RANGE, char *UPLO, int *N, double *A, int *LDA,
             double *VL, double *VU, int *IL, int *IU, double *ABSTOL, int *M,
             double *W, double *Z, int *LDZ, int *ISUPPZ, double *WORK,
             int *LWORK, int *IWORK, int *LIWORK, int *INFO);
}

std::pair<double *, double *> getEigens(const dCSRmat *L, int num_eigens) {
  int n = L->row;

  char jobz = 'V', range = 'I', uplo = 'U';
  int il = 1, m, isuppz[2 * num_eigens], lwork = 26 * n, liwork = 10 * n, info;
  int *iwork = new int[liwork];
  double vl, vu, abstol = 0;
  double *w = new double[n], *z = new double[n * num_eigens],
         *work = new double[lwork];

  double *a = new double[n * n]();
  for (auto i = 0; i < n; ++i) {
    for (auto ind = L->IA[i]; ind < L->IA[i + 1]; ++ind) {
      int j = L->JA[ind];
      if (i <= j) {
        a[i + j * n] = L->val[ind];
      }
    }
  }

  dsyevr_(&jobz, &range, &uplo, &n, a, &n, &vl, &vu, &il, &num_eigens, &abstol,
          &m, w, z, &n, isuppz, work, &lwork, iwork, &liwork, &info);

  assert(m == num_eigens);

  delete[] a;
  delete[] iwork;
  delete[] work;

  return {w, z};
}

std::vector<REAL *> getRandomSmoothVectors(const dCSRmat *L, int num) {
  int n = L->row;
  assert(num <= n);
  int num_eigens = (int)ceil(log2(n));

  auto eigens = getEigens(L, num_eigens);
  delete[] std::get<0>(eigens);
  auto z = std::get<1>(eigens);

  std::default_random_engine generator;
  std::normal_distribution<double> normal(1.0, 1.0);

  std::vector<REAL *> vectors;
  for (int i = 0; i < num; ++i) {
    REAL *v = (REAL *)calloc(n, sizeof(REAL));
    for (int j = 0; j < num_eigens; ++j) {
      double weight = normal(generator);
      for (int k = 0; k < n; ++k)
        v[k] += weight * z[k + j * n];
    }
    vectors.push_back(v);
  }

  delete[] z;

  return vectors;
}

// DIMACS10
// Netherlands_osm

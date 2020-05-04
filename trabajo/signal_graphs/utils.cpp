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

std::vector<REAL *> getRandomSmoothVectors(const dCSRmat *L, int num) {
  int n = L->row;
  assert(num <= n);
  int num_eigen = (int)ceil(log2(n));

  char jobz = 'V', range = 'I', uplo = 'U';
  int il = 1, m, isuppz[2 * num_eigen], lwork = 26 * n, liwork = 10 * n, info;
  int *iwork = new int[liwork];
  double vl, vu, abstol = 1e-1;
  double *w = new double[n], *z = new double[n * num_eigen],
         *work = new double[lwork];

  double *a = new double[n * n]();
  for (auto i = 0; i < n; ++i) {
    for (auto ind = L->IA[i]; ind < L->IA[i + 1]; ++ind) {
      int j = L->JA[ind];
      if (i <= j) {
        a[i * n + j] = L->val[ind];
      }
    }
  }

  dsyevr_(&jobz, &range, &uplo, &n, a, &n, &vl, &vu, &il, &num_eigen, &abstol,
          &m, w, z, &n, isuppz, work, &lwork, iwork, &liwork, &info);
  delete[] a;
  delete[] iwork;
  delete[] w;
  delete[] work;

  std::default_random_engine generator;
  std::normal_distribution<double> normal(1.0, 1.0);

  std::vector<REAL *> vectors;
  for (int i = 0; i < num; ++i) {
    REAL *v = (REAL *)calloc(n, sizeof(REAL));
    for (int j = 0; j < num_eigen; ++j) {
      double weight = normal(generator);
      for (int k = 0; k < n; ++k)
        v[k] += weight * z[k * num_eigen + j];
    }
    vectors.push_back(v);
  }

  delete[] z;

  return vectors;
}

// DIMACS10
// Netherlands_osm

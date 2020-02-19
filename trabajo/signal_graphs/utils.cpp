#include "hazmath_include.h"
#include <algorithm>
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
  char jobz = 'V', range = 'I', uplo = 'U';
  int n = L->row;
  num = std::min(n, num);
  int il = 1, m, isuppz[2 * num], lwork = 26 * n, liwork = 10 * n,
      iwork[liwork], info;
  double vl, vu, abstol = 1e-1, w[num], z[n][num], work[lwork];

  double *a = (double *)calloc(n * n, sizeof(double));
  for (auto i = 0; i < n; ++i) {
    for (auto ind = L->IA[i]; ind < L->IA[i + 1]; ++ind) {
      int j = L->JA[ind];
      if (i <= j) {
        a[i * n + j] = L->val[ind];
      }
    }
  }

  dsyevr_(&jobz, &range, &uplo, &n, a, &n, &vl, &vu, &il, &num, &abstol, &m, w,
          *z, &n, isuppz, work, &lwork, iwork, &liwork, &info);
  free(a);

  std::default_random_engine generator;
  std::normal_distribution<double> normal(1.0, 1.0);

  std::vector<REAL *> vectors;
  for (int i = 0; i < num; ++i) {
    REAL *v = (REAL *)calloc(n, sizeof(REAL));
    for (int j = 0; j < n; ++j) {
      array_axpy(n, normal(generator), z[j], v);
    }
    vectors.push_back(v);
  }
  return vectors;
}

// DIMACS10
// Netherlands_osm

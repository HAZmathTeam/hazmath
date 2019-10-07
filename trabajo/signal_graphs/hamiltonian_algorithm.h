#include <algorithm>
#include <cstdlib>
#include <numeric>
#include <queue>
#include <random>
#include <vector>
#include "hazmath_include.h"

dCSRmat* getRandomProlongation(int n, int N, int seed) {
  assert(n < N);

  std::vector<int> samples(n);
  std::iota(samples.begin(), samples.end(), 0);
  std::shuffle(
      samples.begin(),
      samples.end(),
      std::default_random_engine(seed));
  std::sort(samples.begin(), samples.begin() + N - n);

  dCSRmat *prolong_mat = (dCSRmat*)malloc(sizeof(dCSRmat));
  *prolong_mat = dcsr_create(N, n, N);
  auto sample_it = samples.begin();
  for (int row = 0, vertex = 0; ; ++vertex) {
    prolong_mat->IA[row] = row;
    if (row >= N) {
      break;
    }
    prolong_mat->JA[row] = vertex;
    prolong_mat->val[row] = 1;
    if (vertex == *sample_it) {
      ++row;
      prolong_mat->IA[row] = row;
      prolong_mat->JA[row] = vertex;
      prolong_mat->val[row] = 1;
      ++sample_it;
    }
    ++row;
  }

  return prolong_mat;
}

REAL* prolongate(
    REAL* v,
    const std::vector<int> permutation,
    dCSRmat* prolong_mat) {
  int n = permutation.size();

  REAL permuted_v[n];
  for (int i = 0; i < n; ++i) {
    permuted_v[i] = v[permutation[i]];
  }

  REAL* prolongated_v = (REAL*)calloc(n, sizeof(REAL));
  dcsr_mxv(prolong_mat, permuted_v, prolongated_v);

  return prolongated_v;
}

std::vector<REAL*> getWalshBasis(int L) {
  assert(L >= 0);

  const int N = 1 << L;
  std::vector<REAL*> basis(N);
  basis[0] = (REAL*)calloc(N, sizeof(REAL));
  for (int j = 0; j < N; ++j) {
    basis[0][j] = 1;
  }

  std::queue<std::pair<int, int>> q;
  q.push({1, 0});
  --L;
  while (!q.empty()) {
    int i = q.front().first;
    int l = q.front().second;
    const int b = i - (1 << l);
    basis[i] = (REAL*)calloc(N, sizeof(REAL));
    for (int j = 0; j < N; ++j) {
      basis[i][j] = basis[b][j] * ((j >> (L - l)) % 2 ? -1 : 1);
    }

    if (l < L) {
      i <<= 1;
      ++l;
      q.push({i, l});
      q.push({i + 1, l});
    }
    q.pop();
  }

  return basis;
}

template<typename T>
void deleteArray(const std::vector<T*>& array) {
  for (auto e : array) {
    delete e;
  }
}

REAL* approximate(REAL* v, int N, int k) {

}

REAL*
compress(REAL* v, const std::vector<int> permutation, int k, int seed) {
  int n = permutation.size();
  int N = 1;
  int L = 0;
  while (N < n) {
    N <<= 1;
    ++L;
  }
  dCSRmat* prolong_mat = getRandomProlongation(n, N, seed);
  REAL* prolongated_v = prolongate(v, permutation, prolong_mat);
  REAL* k_term_approximation = approximate(prolongated_v, N, k);
  // REAL* approximation =
  //     project(k_term_approximation, prolong_mat, permutation);

  dcsr_free(prolong_mat);
  delete prolongated_v;
  // delete k_term_approximation;
}

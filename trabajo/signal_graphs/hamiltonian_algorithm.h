#include <algorithm>
#include <cstdlib>
#include <numeric>
#include <random>
#include <vector>
#include "hazmath_include.h"

dCSRmat* getRandomProlongation(int n, int N, int seed) {
  assert(n < N);

  std::vector<int> samples(n);
  std::iota(samples.begin(), samples.end(), 0);
  std::shuffle(
    samples.begin(), samples.end(), std::default_random_engine(seed));
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

double*
compress(double* v, const std::vector<int> permutation, int k, int seed) {
  int n = permutation.size();
  int N = 1;
  while (N < n) {
    N <<= 1;
  }
  dCSRmat* prolong_mat = getRandomProlongation(n, N, seed);
  // double* prolongated_v = prolongate(v, permutation, prolong_mat);
  // double* k_term_approximation = approximate(prolongated_v, N, k);
  // double* approximation =
  //     project(k_term_approximation, prolong_mat, permutation);
}

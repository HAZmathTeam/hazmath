#include "hazmath_include.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <numeric>
#include <queue>
#include <random>
#include <vector>

std::vector<int> getRandomProlongation(int n, int N, int seed) {
  assert(n <= N);

  std::vector<int> samples(N);
  std::iota(samples.begin(), samples.begin() + n, 0);
  std::shuffle(samples.begin(), samples.begin() + n,
               std::default_random_engine(seed));
  std::iota(samples.begin() + N - n, samples.end(), 0);
  std::sort(samples.begin(), samples.end());

  return samples;
}

std::vector<REAL> prolongate(const std::vector<REAL> &v,
                             const std::vector<int> &permutation,
                             const std::vector<int> &samples) {
  int n = permutation.size();
  int N = samples.size();

  std::vector<REAL> permuted_v(n);
  for (int i = 0; i < n; ++i) {
    permuted_v[i] = v[permutation[i]];
  }

  std::vector<REAL> prolongated_v(N);
  for (int i = 0; i < N; ++i) {
    prolongated_v[i] = permuted_v[samples[i]];
  }

  return prolongated_v;
}

std::vector<REAL *> getWalshBasis(int L) {
  assert(L >= 0);

  const int N = 1 << L;
  std::vector<REAL *> basis(N);
  basis[0] = (REAL *)malloc(N * sizeof(REAL));
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
    basis[i] = (REAL *)malloc(N * sizeof(REAL));
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

template <typename T> void deleteArray(const std::vector<T *> &array) {
  for (auto e : array) {
    delete e;
  }
}

std::vector<REAL> approximate(const std::vector<REAL> &v, int L, int k) {
  int N = 1 << L;
  auto basis = getWalshBasis(L);
  std::vector<REAL> inner_products;
  inner_products.reserve(N);
  for (auto b : basis) {
    inner_products.push_back(array_dotprod(N, v.data(), b));
  }

  std::priority_queue<REAL, std::vector<REAL>, std::greater<REAL>> heap;
  for (auto product : inner_products) {
    product = abs(product);
    if (heap.size() < k) {
      heap.push(product);
    } else if (product > heap.top()) {
      heap.pop();
      heap.push(product);
    }
  }
  auto kth_largest = heap.top();

  std::vector<REAL> compression(N, 0.0);
  for (int i = 0; i < N; ++i) {
    auto product = inner_products[i];
    if (abs(product) >= kth_largest) {
      array_axpy(N, product / N, basis[i], compression.data());
    }
  }

  deleteArray(basis);
  return compression;
}

std::vector<REAL> project(const std::vector<REAL> &k_term_approximation,
                          const std::vector<int> &samples,
                          const std::vector<int> &permutation) {
  int n = permutation.size();

  std::vector<REAL> projection;
  projection.reserve(n);
  projection.push_back({k_term_approximation[0]});
  int curr_count = 1;
  for (int i = 1; i < k_term_approximation.size(); ++i) {
    if (samples[i] != samples[i - 1]) {
      projection.back() /= curr_count;
      projection.push_back(k_term_approximation[i]);
      curr_count = 1;
    } else {
      projection.back() += k_term_approximation[i];
      ++curr_count;
    }
  }

  std::vector<REAL> approximation(n);
  for (int i = 0; i < n; ++i) {
    approximation[permutation[i]] = projection[i];
  }

  return approximation;
}

std::vector<REAL> compress(const std::vector<REAL> &v,
                           const std::vector<int> permutation, int k,
                           int seed) {
  int n = permutation.size();
  int N = 1;
  int L = 0;
  while (N < n) {
    N <<= 1;
    ++L;
  }

  const auto samples = getRandomProlongation(n, N, seed);
  const std::vector<REAL> prolongated_v = prolongate(v, permutation, samples);
  const std::vector<REAL> k_term_approximation =
      approximate(prolongated_v, L, k);
  std::vector<REAL> approximation =
      project(k_term_approximation, samples, permutation);

  return approximation;
}

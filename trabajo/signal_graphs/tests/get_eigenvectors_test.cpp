#include "../graph.h"
#include "../utils.h"
#include "utils.h"
#include <algorithm>

std::vector<double> normalize(std::vector<double> &&vector) {
  auto sqr_sum = 0;
  for (auto e : vector) {
    sqr_sum += e * e;
  }
  auto norm = sqrt(sqr_sum);
  for (auto &e : vector) {
    e /= norm;
  }
  return std::move(vector);
}

int main() {
  const Graph graph("graphs/simple.mtx");
  dCSRmat *L = graph.getLaplacian();
  auto n = graph.size();

  auto eigenvalues = std::vector<double>({0, 3 - sqrt(2), 3, 3 + sqrt(2), 5});
  auto eigenvectors = std::vector<std::vector<double>>(
      {normalize({1, 1, 1, 1, 1}),
       normalize({0, -5 - 3 * sqrt(2), -1 - 2 * sqrt(2), 1 + 2 * sqrt(2),
                  (5 + 3 * sqrt(2))}),
       normalize({0, 1, -1, -1, 1}),
       normalize({0, 5 - 3 * sqrt(2), 1 - 2 * sqrt(2), -1 + 2 * sqrt(2),
                  -5 + 3 * sqrt(2)}),
       normalize({-4, 1, 1, 1, 1})});

  auto eigens = getEigens(L, n);
  dcsr_free(L);
  auto w = std::get<0>(eigens);
  auto z = std::get<1>(eigens);

  assertArraysEqual(w, eigenvalues);
  for (auto i = 0; i < n; ++i) {
    try {
      assertArraysEqual(z + i * n, eigenvectors[i], 1e-1);
    } catch (const std::runtime_error &) {
      std::transform(eigenvectors[i].cbegin(), eigenvectors[i].cend(),
                     eigenvectors[i].begin(), std::negate<double>());
      assertArraysEqual(z + i * n, eigenvectors[i], 1e-1);
    }
  }

  delete[] w;
  delete[] z;

  return 0;
}

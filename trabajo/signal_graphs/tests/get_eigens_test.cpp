#include "../graph.h"
#include "../utils.h"
#include "utils.h"
#include <algorithm>

int main() {
  const Graph graph("graphs/simple.mtx");
  dCSRmat *L = graph.getLaplacian();
  auto n = graph.size();

  auto eigenvalues = std::vector<double>({0, 3 - sqrt(2), 3, 3 + sqrt(2), 5});
  // Eigenvectors:
  // {1, 1, 1, 1, 1}
  // {0, -5 - 3 * sqrt(2), -1 - 2 * sqrt(2), 1 + 2 * sqrt(2), 5 + 3 * sqrt(2)}
  // {0, 1, -1, -1, 1}
  // {0, 5 - 3 * sqrt(2), 1 - 2 * sqrt(2), -1 + 2 * sqrt(2), -5 + 3 * sqrt(2)}
  // {-4, 1, 1, 1, 1}

  auto eigens = getEigens(L, n);
  auto w = std::get<0>(eigens);

  dcsr_free(L);
  delete[] std::get<1>(eigens);

  assertArraysEqual(w, eigenvalues);

  delete[] w;

  return 0;
}

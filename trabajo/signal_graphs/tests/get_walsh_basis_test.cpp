#include "../hamiltonian_algorithm.h"
#include "utils.h"

int main() {
  const std::vector<std::vector<std::vector<double>>> expected_bases(
      {{{1}},
       {{1, 1}, {1, -1}},
       {{1, 1, 1, 1}, {1, 1, -1, -1}, {1, -1, 1, -1}, {1, -1, -1, 1}},
       {{1, 1, 1, 1, 1, 1, 1, 1},
        {1, 1, 1, 1, -1, -1, -1, -1},
        {1, 1, -1, -1, 1, 1, -1, -1},
        {1, 1, -1, -1, -1, -1, 1, 1},
        {1, -1, 1, -1, 1, -1, 1, -1},
        {1, -1, 1, -1, -1, 1, -1, 1},
        {1, -1, -1, 1, 1, -1, -1, 1},
        {1, -1, -1, 1, -1, 1, 1, -1}}});
  for (int l = 0; l < expected_bases.size(); ++l) {
    auto basis = getWalshBasis(l);
    auto &expected_basis = expected_bases[l];
    for (int i = 0; i < expected_basis.size(); ++i) {
      assertArraysEqual(basis[i], expected_basis[i]);
    }
    deleteArray(basis);
  }

  return 0;
}
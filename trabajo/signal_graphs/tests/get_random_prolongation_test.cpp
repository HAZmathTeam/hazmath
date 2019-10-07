#include "utils.h"
#include "../hamiltonian_algorithm.h"

int main(int argc, char *argv[]) {
  auto* prolong_mat = getRandomProlongation(3, 5, 0);
  assertArraysEqual(prolong_mat->IA, std::vector<int>({0, 1, 2, 3, 4, 5}));
  assertArraysEqual(prolong_mat->JA, std::vector<int>({0, 0, 1, 2, 2}));
  assertArraysEqual(prolong_mat->val, std::vector<double>({1, 1, 1, 1, 1}));

  dcsr_free(prolong_mat);
  return 0;
}

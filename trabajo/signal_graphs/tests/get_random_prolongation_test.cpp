#include "../hamiltonian_algorithm.h"

int main(int argc, char *argv[]) {
  auto* prolong_mat = getRandomProlongation(3, 5, 0);

  std::vector<int> IA({0, 1, 2, 3, 4, 5});
  std::vector<int> JA({0, 0, 1, 2, 2});
  std::vector<int> val({1, 1, 1, 1, 1});
  for (int i = 0; i < IA.size(); ++i) {
    assert(prolong_mat->IA[i] == IA[i]);
  }
  for (int i = 0; i < JA.size(); ++i) {
    assert(prolong_mat->JA[i] == JA[i]);
    assert(prolong_mat->val[i] == val[i]);
  }

  dcsr_free(prolong_mat);

  return 0;
}

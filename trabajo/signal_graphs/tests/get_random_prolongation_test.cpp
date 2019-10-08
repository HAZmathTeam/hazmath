#include "utils.h"
#include "../hamiltonian_algorithm.h"

int main(int argc, char *argv[]) {
  auto samples = getRandomProlongation(3, 5, 0);
  assertArraysEqual(samples.data(), std::vector<int>({0, 0, 1, 2, 2}));

  return 0;
}

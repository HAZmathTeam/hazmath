#include "../graph.h"
#include "../algorithm.h"

int main(int argc, char *argv[]) {
  auto&& path1 = get_hamiltonian_path(Tree(1, {Tree(2, {Tree(3)})}));
  assert(path1 == std::vector<int>({1, 3, 2}));

  auto&& path2 = get_hamiltonian_path(Tree(1, {Tree(2), Tree(3)}));
  assert(
      path2 == std::vector<int>({1, 2, 3}) ||
      path2 == std::vector<int>({1, 3, 2}));

  return 0;
}

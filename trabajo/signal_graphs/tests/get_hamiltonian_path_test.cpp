#include "../graph.h"
#include "../algorithm.h"

int main(int argc, char *argv[]) {
  auto&& path1 =
    get_hamiltonian_path(new Tree(1, {new Tree(2, {new Tree(3)})}));
  assert(path1 == std::vector<int>({1, 3, 2}));

  auto&& path2 = get_hamiltonian_path(new Tree(1, {new Tree(2), new Tree(3)}));
  assert(
      path2 == std::vector<int>({1, 2, 3}) ||
      path2 == std::vector<int>({1, 3, 2}));

  Graph graph("graphs/simple.mtx");
  auto&& path3 = graph.GetHamiltonianPath(1);
  assert(path3 == std::vector<int>({2, 4, 0, 3, 1}));

  return 0;
}

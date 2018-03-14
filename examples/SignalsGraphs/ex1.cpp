/* Example to apply data compression algorithm on graphs
 * Usage:
 *   ./ex1 graphs/power.mtx 1000
 */
#include <string>
#include "graph.hpp"
#include "algorithm.hpp"

using namespace std;

int main(int argc, char *argv[]) {
  assert(argc > 1);
  int threshold = 100;
  if (argc > 2) {
    threshold = stoi(argv[2]);
  }
  double p = 1.0;
  if (argc > 3) {
    p = stod(argv[3]);
  }

  dCSRmat *A;
  vector<dCSRmat *> Qj_array;
  vector<int> Nj_array;
  setup_hierarchy(argv[1], A, Qj_array, Nj_array);
  int n = A->row;
  if (threshold > n - 1) {
    threshold = n - 1;
  }
  REAL *v2 = (REAL *)malloc(sizeof(REAL)*n);
  REAL *v3 = (REAL *)malloc(sizeof(REAL)*n);
  comp_decomp(NULL, A, Qj_array, Nj_array, threshold, p, v2, v3);

  dcsr_free(A);
  for (auto Qj : Qj_array) {
    dcsr_free(Qj);
  }
  free(v2);
  free(v3);
  return 0;
}

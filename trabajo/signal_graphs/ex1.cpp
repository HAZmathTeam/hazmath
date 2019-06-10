/* Example to apply data compression algorithm on graphs
 * Usage:
 *   ./ex1 graphs/power.mtx -k 100 -p 1.0
 */
#include <iostream>
#include <string>
#include "graph.hpp"
#include "algorithm.hpp"

using namespace std;

int main(int argc, char *argv[]) {
  int threshold = 100;
  double p = 1.0;
  int c;
  opterr = 0;

  while ((c = getopt(argc, argv, "k:p:")) != -1) {
    switch (c) {
      case 'k':
        threshold = stoi(optarg);
        break;
      case 'p':
        p = stod(optarg);
        break;
      case '?':
        if (optopt == 'k' || optopt == 'p') {
          fprintf(stderr, "Option -%c requires an argument.\n", optopt);
        }
        else if (isprint(optopt)) {
          fprintf(stderr, "Unknown option `-%c'.\n", optopt);
        }
        else {
          fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
        }
        return 1;
      default:
        abort();
    }
  }
  if (optind != argc - 1) {
    cerr << "Too many arguments." << endl;
    return 1;
  }

  dCSRmat *A;
  vector<dCSRmat *> Qj_array;
  vector<int> Nj_array;
  setup_hierarchy(argv[optind], A, Qj_array, Nj_array);
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

/* Example to apply data compression algorithm on graphs
 * Usage:
 *   ./ex1 -k 100 -p 1.0 graphs/power.mtx
 */
#include "algorithm.h"
#include "graph.h"
#include "utils.h"
#include <iostream>
#include <string>

using namespace std;

int main(int argc, char *argv[]) {
  int largestK = 100;
  double p = 1.0;
  int c;
  opterr = 0;

  while ((c = getopt(argc, argv, "k:p:")) != -1) {
    switch (c) {
    case 'k':
      largestK = stoi(optarg);
      break;
    case 'p':
      p = stod(optarg);
      break;
    case '?':
      if (optopt == 'k' || optopt == 'p') {
        fprintf(stderr, "Option -%c requires an argument.\n", optopt);
      } else if (isprint(optopt)) {
        fprintf(stderr, "Unknown option `-%c'.\n", optopt);
      } else {
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

  const Graph graph(argv[optind]);

  dCSRmat *L = graph.getLaplacian();
  std::vector<REAL *> vectors = getRandomSmoothVectors(L);
  dcsr_free(L);

  if (largestK > graph.size() - 1) {
    largestK = graph.size() - 1;
  }
  std::cout << std::endl << std::endl;
  ConnectionMatchingWalsh().compAndDecomp(graph, vectors, largestK, p);
  std::cout << std::endl << std::endl;
  DegreeMatchingWalsh().compAndDecomp(graph, vectors, largestK, p);
  std::cout << std::endl << std::endl;
  ConnectionMatchingGtbwt().compAndDecomp(graph, vectors, largestK, p);
  std::cout << std::endl << std::endl;
  DegreeMatchingGtbwt().compAndDecomp(graph, vectors, largestK, p);
  std::cout << std::endl << std::endl;
  HamiltonianAlgorithm().compAndDecomp(graph, vectors, largestK);

  for (REAL *v : vectors) {
    free(v);
  }
  return 0;
}

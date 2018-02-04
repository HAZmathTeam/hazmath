/* Example to apply data compression algorithm on graphs
 * Usage:
 *   ./ex1 graphs/power.mtx 1000
 */
#include <iostream>
#include <numeric>
#include <algorithm>
#include <functional>
#include "graph.hpp"

using namespace std;

extern "C" {
  void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda, double *w,
              double *work, int *lwork, int *info);
}

void setup_hierarchy(const char *file, dCSRmat *&A, vector<dCSRmat *> &Qj_array,
    vector<int> &Nj_array);

void comp_decomp(int argc, char *argv[], double *v, dCSRmat *A,
    const vector<dCSRmat *> &Qj_array, const vector<int> &Nj_array,
    double* v2, double *v3);

int main(int argc, char *argv[]) {
  assert(argc > 1);

  dCSRmat *A;
  vector<dCSRmat *> Qj_array;
  vector<int> Nj_array;
  setup_hierarchy(argv[1], A, Qj_array, Nj_array);
  int n = A->row;
  REAL *v2 = (REAL *)malloc(sizeof(REAL)*n);
  REAL *v3 = (REAL *)malloc(sizeof(REAL)*n);
  comp_decomp(argc-1, argv+1, NULL, A, Qj_array, Nj_array, v2, v3);

  dcsr_free(A);
  for (auto Qj : Qj_array) {
    dcsr_free(Qj);
  }
  free(v2);
  free(v3);
  return 0;
}

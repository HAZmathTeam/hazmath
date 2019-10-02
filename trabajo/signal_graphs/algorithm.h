#ifndef SIGNALS_GRAPH_ALGORITHM
#define SIGNALS_GRAPH_ALGORITHM

#ifndef EXTERN_C
#define EXTERN_C
extern "C" {
  #include "hazmath.h"
}
#endif

#include "graph.h"

void setup_hierarchy(const char *file, dCSRmat *&A,
    std::vector<dCSRmat *> &Qj_array, std::vector<int> &Nj_array);

void comp_decomp(double *v, dCSRmat *A, const std::vector<dCSRmat *> &Qj_array,
    const std::vector<int> &Nj_array, int threshold, double p, double* v2,
    double *v3);

std::vector<int> get_hamiltonian_path(Tree&& tree);

#endif

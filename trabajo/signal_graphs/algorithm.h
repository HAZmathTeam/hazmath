#ifndef SIGNALS_GRAPH_ALGORITHM
#define SIGNALS_GRAPH_ALGORITHM

#include "graph.h"
#include "hazmath_include.h"

void setupHierarchy(const char *file, dCSRmat *&A,
                    std::vector<dCSRmat *> &Qj_array,
                    std::vector<int> &Nj_array);

void compAndDecomp(double *v, dCSRmat *A,
                   const std::vector<dCSRmat *> &Qj_array,
                   const std::vector<int> &Nj_array, int threshold, double p,
                   double *v2, double *v3);

std::vector<int> getHamiltonianPath(Tree *tree);

#endif

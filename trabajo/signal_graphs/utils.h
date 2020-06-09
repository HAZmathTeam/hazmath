#ifndef SIGNALS_GRAPH_UTILS
#define SIGNALS_GRAPH_UTILS

#include "hazmath_include.h"
#include <utility>

REAL *initializeRhs(dCSRmat *A, int num_iterations = 100);

std::pair<double *, double *> getEigens(const dCSRmat *L, int num_eigens);

std::vector<REAL *> getRandomSmoothVectors(const dCSRmat *L, int num = 20);

#endif

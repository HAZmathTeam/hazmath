#ifndef SIGNALS_GRAPH_UTILS
#define SIGNALS_GRAPH_UTILS

#include "hazmath_include.h"

REAL *initializeRhs(dCSRmat *A, int num_iterations = 100);

std::vector<REAL *> getRandomSmoothVectors(const dCSRmat *L, int num = 20);

#endif

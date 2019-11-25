#ifndef SIGNALS_GRAPH_HAMILTONIAN
#define SIGNALS_GRAPH_HAMILTONIAN

#include "graph.h"
#include "hazmath_include.h"

std::vector<int> getHamiltonianPath(Tree *tree);

std::vector<int> getRandomProlongation(int n, int N, int seed);

std::vector<REAL> prolongate(const std::vector<REAL> &v,
                             const std::vector<int> &permutation,
                             const std::vector<int> &samples);

std::vector<REAL *> getWalshBasis(int L);

template <typename T> void deleteArray(const std::vector<T *> &array);

std::vector<REAL> approximate(const std::vector<REAL> &v, int L, int k);

std::vector<REAL> project(const std::vector<REAL> &k_term_approximation,
                          const std::vector<int> &samples,
                          const std::vector<int> &permutation);

std::vector<REAL> compress(const std::vector<REAL> &v,
                           const std::vector<int> permutation, int k, int seed);

#endif

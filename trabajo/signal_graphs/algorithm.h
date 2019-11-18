#ifndef SIGNALS_GRAPH_ALGORITHM
#define SIGNALS_GRAPH_ALGORITHM

#include "graph.h"
#include "hazmath_include.h"

class Algorithm {
public:
  virtual int numBlocks(int numBlocks) const = 0;
};

class Adaptive : public Algorithm {
public:
  int numBlocks(int numBlocks) const { return numBlocks; }
};

class Gtbwt : public Algorithm {
public:
  int numBlocks(int numBlocks) const { return 1; }
};

void setupHierarchy(Graph graph, std::vector<dCSRmat *> &Qj_array,
                    std::vector<int> &Nj_array, const Algorithm &algorithm);

void compAndDecomp(int n, double *v, const std::vector<dCSRmat *> &Qj_array,
                   int largestK, double p, double *v2);

void compAndDecomp(int n, double *v, const std::vector<dCSRmat *> &Qj_array,
                   const std::vector<int> &Nj_array, int largestK, double p,
                   double *v3);

std::vector<int> getHamiltonianPath(Tree *tree);

void compAndDecomp(const Graph &graph, REAL *v, const int largestK,
                   const double p, const Algorithm &algorithm);

#endif

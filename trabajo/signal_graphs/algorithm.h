#ifndef SIGNALS_GRAPH_ALGORITHM
#define SIGNALS_GRAPH_ALGORITHM

#include "graph.h"
#include "hazmath_include.h"

class Algorithm {
public:
  virtual bool isAdaptive() const { return false; }

  void setupHierarchy(Graph graph, std::vector<dCSRmat *> &Qj_array,
                      std::vector<int> &Nj_array) const;

  void compAndDecomp(int n, double *v, const std::vector<dCSRmat *> &Qj_array,
                     int largestK, double p, double *v2) const;

  void compAndDecomp(int n, double *v, const std::vector<dCSRmat *> &Qj_array,
                     const std::vector<int> &Nj_array, int largestK, double p,
                     double *v3) const;

  void compAndDecomp(const Graph &graph, REAL *v, const int largestK,
                     const double p) const;

private:
  virtual int getNumBlocks(int numBlocks) const = 0;
};

class Adaptive : public Algorithm {
public:
  bool isAdaptive() const { return true; }

private:
  int getNumBlocks(int numBlocks) const { return numBlocks; }
};

class Gtbwt : public Algorithm {
private:
  int getNumBlocks(int numBlocks) const { return 1; }
};

std::vector<int> getHamiltonianPath(Tree *tree);

#endif

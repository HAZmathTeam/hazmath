#ifndef SIGNALS_GRAPH_ALGORITHM
#define SIGNALS_GRAPH_ALGORITHM

#include "graph.h"
#include "hazmath_include.h"
#include <string>

class Algorithm {
public:
  virtual void compAndDecomp(const Graph &graph, REAL *v, const int largestK,
                             const double p) const = 0;
};

class AggregationBasedAlgorithm : public Algorithm {
public:
  virtual bool isAdaptive() const { return false; }

  void setupHierarchy(Graph graph, std::vector<dCSRmat *> &Qj_array,
                      std::vector<int> &Nj_array, bool weighted = true) const;

  void compAndDecomp(int n, double *v, const std::vector<dCSRmat *> &Qj_array,
                     int largestK, double *v2) const;

  void compAndDecomp(int n, double *v, const std::vector<dCSRmat *> &Qj_array,
                     const std::vector<int> &Nj_array, int largestK, double p,
                     double *v3) const;

  void compAndDecomp(const Graph &graph, REAL *v, const int largestK,
                     const double p) const;

private:
  virtual std::string matchingAlgorithm() const = 0;

  virtual std::string compressionAlgorithm() const = 0;

  virtual void doMatching(Graph &graph, Graph *c_graph) const = 0;

  virtual int getNumBlocks(int numBlocks) const = 0;
};

class ConnectionBasedMatching : virtual public AggregationBasedAlgorithm {
private:
  std::string matchingAlgorithm() const { return "connection-based"; };

  void doMatching(Graph &graph, Graph *c_graph) const {
    graph.doConnectionBasedMatching(c_graph);
  }
};

class DegreeBasedMatching : virtual public AggregationBasedAlgorithm {
private:
  std::string matchingAlgorithm() const { return "degree-based"; };

  void doMatching(Graph &graph, Graph *c_graph) const {
    graph.doDegreeBasedMatching(c_graph);
  }
};

class Walsh : virtual public AggregationBasedAlgorithm {
public:
  bool isAdaptive() const { return true; }

private:
  std::string compressionAlgorithm() const { return "adaptive"; };

  int getNumBlocks(int numBlocks) const { return numBlocks; }
};

class Gtbwt : virtual public AggregationBasedAlgorithm {
private:
  std::string compressionAlgorithm() const { return "GTBWT"; };

  int getNumBlocks(int numBlocks) const { return 1; }
};

class ConnectionMatchingWalsh : public ConnectionBasedMatching, Walsh {};

class DegreeMatchingWalsh : public DegreeBasedMatching, Walsh {};

class ConnectionMatchingGtbwt : public ConnectionBasedMatching, Gtbwt {};

class DegreeMatchingGtbwt : public DegreeBasedMatching, Gtbwt {};

class HamiltonianAlgorithm {
public:
  void compAndDecomp(const Graph &graph, REAL *v, const int largestK) const;
};

#endif

#ifndef SIGNALS_GRAPH
#define SIGNALS_GRAPH

#include <vector>

extern "C" {
  #include "hazmath.h"
}

class Graph {
private:
  // Adjacecy matrix of vertices
  iCSRmat *A;

  // The aggregates after matching
  std::vector<std::vector<int>> aggregates;

public:
  Graph(): A(NULL) {}

  Graph(const char* filename);

  // Number of vertices
  int Size() const { return A->row; }

  // Perform matching algorithm and construct the coarse graph
  void DoMatching(Graph *c_graph);

  // Get number of aggregates in the graph
  int NumOfAggregates() const { return aggregates.size(); }

  // Get the vertices in an aggregate
  void GetAggregate(int i, std::vector<int> *vertices) const;

  // Get the number of vertices in an aggregate
  int GetAggregateSize(int i) const { return aggregates[i].size(); }

  // Get the graph Laplacian
  dCSRmat *GetWeightedLaplacian() const;
};

#endif

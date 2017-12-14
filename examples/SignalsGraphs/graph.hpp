#ifndef SIGNALS_GRAPH
#define SIGNALS_GRAPH

#include <vector>
#include <utility>

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

  Graph(const Graph &other) {
    A = (iCSRmat *)malloc(sizeof(iCSRmat));
    *A = icsr_create(other.A->row, other.A->col, other.A->nnz);
    icsr_cp(other.A, A);
    aggregates = other.aggregates;
  }

  ~Graph() { icsr_free(A); }

  Graph & operator= (Graph other) {
    std::swap(A, other.A);
    std::swap(aggregates, other.aggregates);

    return *this;
  }

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

#ifndef SIGNALS_GRAPH
#define SIGNALS_GRAPH

#ifndef EXTERN_C
#define EXTERN_C
extern "C" {
  #include "hazmath.h"
}
#endif

#include <vector>

class Graph {
private:
  // Adjacecy matrix of vertices
  iCSRmat *A;

  // The aggregates after matching
  std::vector<std::vector<int>> aggregates;

  std::vector<int> GetNeighbors(int i) const;

public:
  Graph(): A(NULL) {}

  Graph(const char* filename);

  Graph(const Graph& other);

  ~Graph() { icsr_free(A); }

  Graph & operator= (Graph other);

  // Number of vertices
  int Size() const {
    return A->row;
  }

  // Perform matching algorithm and construct the coarse graph
  void DoMatching(Graph *c_graph);

  // Get number of aggregates in the graph
  int NumOfAggregates() const {
    return aggregates.size();
  }

  // Get the vertices in an aggregate
  void GetAggregate(int i, std::vector<int> *vertices) const;

  // Get the number of vertices in an aggregate
  int GetAggregateSize(int i) const {
    return aggregates[i].size();
  }

  // Get the graph Laplacian
  dCSRmat *GetWeightedLaplacian() const;

  // Get Hamiltonian path
  std::vector<int> GetHamiltonianPath(int seed = 0) const;
};

class Tree {
public:
  int vertex;
  std::vector<Tree> children;

  Tree(int vertex, std::vector<Tree>&& children = {})
    : vertex(vertex), children(children) {}
};

#endif

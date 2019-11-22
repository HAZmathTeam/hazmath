#ifndef SIGNALS_GRAPH
#define SIGNALS_GRAPH

#include "hazmath_include.h"
#include <vector>

class Tree {
public:
  int vertex;
  std::vector<Tree *> children;

  Tree(int vertex, std::vector<Tree *> &&children = {})
      : vertex(vertex), children(children) {}
};

class Graph {
private:
  // Adjacency matrix of vertices
  iCSRmat *A;

  // The aggregates after matching
  std::vector<std::vector<int>> aggregates;

  std::vector<int> getNeighbors(int i) const;

public:
  Graph() : A(NULL) {}

  Graph(const char *filename);

  Graph(const Graph &other);

  ~Graph() { icsr_free(A); }

  Graph &operator=(Graph other);

  const iCSRmat *getAdjacencyMat() const { return A; }

  const std::vector<std::vector<int>> &getAggregates() const {
    return aggregates;
  }

  // Number of vertices
  int size() const { return A->row; }

  // Perform matching algorithm and construct the coarse graph
  void doMatching(Graph *c_graph);

  void doDegreeBasedMatching(Graph *c_graph, int seed = 0);

  // Get number of aggregates in the graph
  int numOfAggregates() const { return aggregates.size(); }

  // Get the vertices in an aggregate
  void getAggregate(int i, std::vector<int> *vertices) const;

  // Get the number of vertices in an aggregate
  int getAggregateSize(int i) const { return aggregates[i].size(); }

  // Get the graph Laplacian
  dCSRmat *getWeightedLaplacian() const;

  // Get Hamiltonian path
  std::vector<int> getHamiltonianPath(int seed = 0) const;
};

#endif

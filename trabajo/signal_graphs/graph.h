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
  // Unweighted (weights = 1) adjacency matrix of the graph.
  iCSRmat *A;

  // The aggregates (represented as vertex indices in the current graph) in the
  // subgraph after matching is performed on the current graph.
  std::vector<std::vector<int>> aggregates;

  // Number of vertices in each aggregate from the finest graph.
  std::vector<int> abs_card_;

  std::vector<int> getNeighbors(int i) const;

public:
  Graph() : A(NULL) {}

  Graph(const char *filename);

  Graph(const Graph &other);

  ~Graph() { icsr_free(A); }

  Graph &operator=(Graph &&other);

  const iCSRmat *getAdjacencyMat() const { return A; }

  const std::vector<std::vector<int>> &getAggregates() const {
    return aggregates;
  }

  // Number of vertices
  int size() const { return A->row; }

  // Perform matching algorithm and construct the coarse graph
  void doConnectionBasedMatching(Graph *c_graph);

  void doDegreeBasedMatching(Graph *c_graph, int seed = 0);

  // Get number of aggregates in the graph
  int numOfAggregates() const { return aggregates.size(); }

  // Get the vertices in an aggregate
  void getAggregate(int i, std::vector<int> *vertices) const;

  // Get the number of vertices in an aggregate
  int getAggregateSize(int i) const { return aggregates[i].size(); }

  int getAbsCard(int i) const { return abs_card_[i]; }

  // Get the unweighted graph Laplacian matrix.
  dCSRmat *getLaplacian() const;

  // Get Hamiltonian path
  std::vector<int> getHamiltonianPath(int seed = 0) const;
};

#endif

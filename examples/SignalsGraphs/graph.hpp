#ifndef ADAPTIVE_GRAPH
#define ADAPTIVE_GRAPH

#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <boost/functional/hash.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/mem_fun.hpp>

extern "C" {
#include "hazmath.h"
}

typedef unsigned int Index;
typedef Index AggregateIndex;
typedef Index VertexIndex;
typedef Index EdgeIndex;

class Aggregate;
typedef std::shared_ptr<Aggregate> AggregatePtr;

template <typename T>
size_t UnorderedTwoHash(T a, T b) {
  size_t seed = 0;
  if (a > b) std::swap(a, b);
  boost::hash_combine(seed, a);
  boost::hash_combine(seed, b);
  return seed;
}

struct Edge {
  VertexIndex head;
  VertexIndex tail;

  Edge(VertexIndex head, VertexIndex tail): head(head), tail(tail) {}

  bool operator == (const Edge& other) const {
    return (head == other.head && tail == other.tail) ||
           (head == other.tail && tail == other.head);
  }
};

namespace std {
  template <>
  struct hash<Edge> {
    size_t operator() (const Edge &e) const {
      return UnorderedTwoHash(e.head, e.tail);
    }
  };
}

typedef boost::multi_index_container<
  Edge,
  boost::multi_index::indexed_by<
    boost::multi_index::hashed_unique<
      boost::multi_index::identity<Edge>, std::hash<Edge>
    >,
    boost::multi_index::hashed_non_unique<
      boost::multi_index::member<Edge, VertexIndex, &Edge::head>
    >
  >
> EdgeSet;

struct EdgeVector {
  // using int instead of EdgeIndex for compatibility with mfem::Array<int>
  std::vector<int> edge_indices;
  std::vector<double> values;
};

class Aggregate {
  friend class Graph;

private:
  std::unordered_set<VertexIndex> vertices;
  std::vector<AggregatePtr> sub_aggregates;
  EdgeSet interface_edges;
  std::unordered_map<VertexIndex, EdgeVector> PiHColumns;

public:
  // Number of vertices
  int Size() const { return vertices.size(); }

  bool HasVertex(VertexIndex v) const { return vertices.count(v); }

  // Add a vertex
  inline void AddVertex(VertexIndex v) { vertices.insert(v); }

  // Add an interface edge
  void AddInterfaceEdge(const Edge &e) { interface_edges.insert(e); }

  // Merge a sub aggregate into the current one
  void Merge(const AggregatePtr &agg);

  // Merge two aggregates and return a new one
  friend AggregatePtr Merge(const AggregatePtr &a1, const AggregatePtr &a2);

  // Determine if PiH has been generated
  bool PiHGenerated() const { return !PiHColumns.empty();}
};

struct Interface {
  AggregateIndex head_index;
  AggregateIndex tail_index;
  std::unordered_set<Edge> edges;

  Interface(AggregateIndex h, AggregateIndex t): head_index(h), tail_index(t) {
    if (head_index > tail_index) std::swap(head_index, tail_index);
  }

  // Number of edges in the interface
  int Size() const { return edges.size(); }

  // Hash function for fast lookup in the class InterfaceSet
  size_t hash() const { return UnorderedTwoHash(head_index, tail_index); }

  // Merge a set of edges into the current interface
  void Merge(const std::unordered_set<Edge> &edges) {
    this->edges.insert(edges.begin(), edges.end());
  }

  // Add an edge to the interface
  void Add(const Edge &edge) { edges.insert(edge); }
};

typedef boost::multi_index_container<
  Interface,
  boost::multi_index::indexed_by<
    boost::multi_index::hashed_non_unique<
      boost::multi_index::member<
        Interface, AggregateIndex, &Interface::head_index
      >
    >,
    boost::multi_index::hashed_non_unique<
      boost::multi_index::member<
        Interface, AggregateIndex, &Interface::tail_index
      >
    >,
    boost::multi_index::hashed_unique<
      boost::multi_index::const_mem_fun<Interface, size_t, &Interface::hash>
    >
  >
> BaseInterfaceSet;

class InterfaceSet: public BaseInterfaceSet {
private:
  typedef BaseInterfaceSet::nth_index<2>::type LookupType;
  LookupType &lookup = BaseInterfaceSet::get<2>();

public:
  InterfaceSet() {}

  InterfaceSet(const InterfaceSet &other):
      BaseInterfaceSet(other), lookup(BaseInterfaceSet::get<2>()) {}

  InterfaceSet &operator=(InterfaceSet other) {
    BaseInterfaceSet::swap(other);
    return *this;
  }

  // Look for interface, if not exist then create it
  LookupType::iterator query(AggregateIndex a, AggregateIndex b) {
    auto it = lookup.find(UnorderedTwoHash(a, b));
    return it == lookup.end() ? lookup.insert(Interface(a, b)).first : it;
  }

  // Merge interface into the one pointed by iterator
  void Merge(AggregateIndex a, AggregateIndex b,
             const std::unordered_set<Edge> &edges) {
    LookupType::iterator it = query(a, b);
    Interface interface = *it;
    interface.Merge(edges);
    lookup.replace(it, interface);
  }

  // Add edge to interface pointed by iterator
  void Add(AggregateIndex a, AggregateIndex b, const Edge &edge) {
    LookupType::iterator it = query(a, b);
    Interface interface = *it;
    interface.Add(edge);
    lookup.replace(it, interface);
  }

  // Remove interfaces that contains the aggregate index i
  void remove(AggregateIndex i) {
    remove<0>(i);
    remove<1>(i);
  }

  // N is the index of the hashed index
  template <int N>
  void remove(AggregateIndex i) {
    auto &search = BaseInterfaceSet::get<N>();
    auto it_pair = search.equal_range(i);
    search.erase(it_pair.first, it_pair.second);
  }
};

// The graph class that can represent aggreation hierarchies
class Graph {
private:
  // Number of vertices
  int size;

  // Adjacecy matrix of vertices
  std::vector<std::vector<VertexIndex> > adjacency_table;

  // Vertex to aggregate index map
  std::vector<AggregateIndex> vertex_to_aggregate;

  // Edges of the graph
  std::vector<Edge> edges;

  // Hash table of edges to its index
  std::unordered_map<Edge, EdgeIndex> edge_to_index;

  // Aggregates, stored as shared pointers
  std::vector<AggregatePtr> aggregates;

  // Interfaces
  InterfaceSet interfaces;

  // Bool indicating if reshaping has been performed on the graph
  bool reshaped = false;

  // Generate edges only from an unmatched graph
  void GenerateEdges();

  // Generate the lookup array vertex_to_aggregate
  void GenerateVertexMap();

  // Update the lookup array vertex_to_aggregate from given indices
  void UpdateVertexMap(const std::vector<AggregateIndex> &indices);

  // Generate interfaces for aggregate i
  void GenerateInterfaces(AggregateIndex i);

  // // Generate the columns of the interpolcation Pi_H for an aggregate, the class
  // //  Graph is declared as a friend of the class Aggregate so it has access to
  // //  its private and protected members
  // void GeneratePiHColumns(const AggregatePtr &agg);
  //
  // // The recursive function (implementing dfs) used in constructing the spanning
  // //  trees
  // int TreeSize(VertexIndex root, const AggregatePtr &agg,
  //              std::unordered_set<VertexIndex> &visited, EdgeVector &v);

public:
  Graph(): size(0) {}

  Graph(const char* filename);

  // // Assignment operator
  // Graph &operator=(const Graph &other);

  // Number of vertices
  size_t Size() const { return size; }

  // Perform matching algorithm on the graph
  void DoMatching();

  // Perform matching algorithm and construct the coarse graph
  void DoMatching(Graph *c_graph);

  // Get number of aggregates in the graph
  size_t NumOfAggregates() const { return aggregates.size(); }

  // Get the vertices in an aggregate
  void GetAggregate(AggregateIndex i, std::vector<int> *vertices) const;

  // Get the graph Laplacian
  dCSRmat *GetWeightedLaplacian() const;

  // // Get the coarse grid right-hand side fH = Prol' * f
  // mfem::Vector GetCoarseRHS(const mfem::Vector &f) const;
  //
  // // Prolongate a corase vector to fine vector
  // mfem::Vector Prolongate(const mfem::Vector &uh_c) const;
  //
  // // Get the Phi matrix
  // mfem::SparseMatrix *GetCoarseEdgeSpaceBasis();
  //
  // // Get the localized error estimates
  // mfem::Vector GetLocalEstimates(mfem::Vector f, const mfem::Vector &uh,
  //                                mfem::Vector y, double CP) const;
};

#endif

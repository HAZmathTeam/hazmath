#include <iostream>
#include <fstream>
#include <unordered_map>
#include <cassert>

#include "graph.hpp"
#include "boost/unordered_map.hpp"

using namespace std;

void Aggregate::Merge(const AggregatePtr& agg) {
  this->vertices.insert(agg->vertices.begin(), agg->vertices.end());
  sub_aggregates.push_back(agg);
  // Merge the interface edge sets
  for (auto edge : agg->interface_edges) {
    if (this->interface_edges.get<0>().count(edge)) {
      this->interface_edges.erase(edge);
    }
    else this->interface_edges.insert(edge);
  }
}

AggregatePtr Merge(const AggregatePtr& agg1, const AggregatePtr& agg2) {
  AggregatePtr agg = make_shared<Aggregate>();
  agg->vertices.insert(agg1->vertices.begin(), agg1->vertices.end());
  agg->vertices.insert(agg2->vertices.begin(), agg2->vertices.end());
  agg->sub_aggregates.push_back(agg1);
  agg->sub_aggregates.push_back(agg2);
  // Merge the interface edge sets
  for (auto edge : agg1->interface_edges)
    if (!agg2->interface_edges.get<0>().count(edge))
      agg->interface_edges.insert(edge);
  for (auto edge : agg2->interface_edges)
    if (!agg1->interface_edges.get<0>().count(edge))
      agg->interface_edges.insert(edge);
  return agg;
}

/* Assumption on the mtx file representing a matrix:
 * 1. Vertex numbers start from 1 instead of 0
 * 2. For each edge (i,j) (i < j), the entry (j,i) is not recorded
 * 3. Neighbors of a vertex appear in ascending order
 */
Graph::Graph(const char* filename) {
  ifstream file(filename);
  string line;
  while (getline(file, line) && (line.empty() || line[0] == '%'));
  unsigned nrows, ncols, nnz;
  istringstream iss(line);
  if (!(iss >> nrows >> ncols >> nnz)) {
    throw runtime_error("Error reading matrix sizes! Exiting...");
  }
  if (nrows != ncols) {
    throw logic_error("Not square matrix! Exiting...");
  }
  for (VertexIndex i = 0; i < nrows; ++i) {
    AggregatePtr agg = make_shared<Aggregate>();
    agg->AddVertex(i);
    aggregates.push_back(agg);
  }
  adjacency_table.resize(nrows);
  unsigned count = 0, i, j;
  while (getline(file, line)) {
    ++count;
    istringstream iss(line);
    if (!(iss >> i >> j)) {
      throw runtime_error("Error reading matrix entry! Exiting...");
    }
    if (i == j) continue;
    --i; --j;
    assert(adjacency_table[i].empty() || adjacency_table[i].back() < j);
    adjacency_table[i].push_back(VertexIndex(j));
    aggregates[i]->AddInterfaceEdge(Edge(i, j));
    assert(adjacency_table[j].empty() || adjacency_table[j].back() < i);
    adjacency_table[j].push_back(VertexIndex(i));
    aggregates[j]->AddInterfaceEdge(Edge(j, i));

    interfaces.Add(i, j, Edge(i, j));
  }
  assert(count == nnz);

  size = nrows;
  vertex_to_aggregate.resize(size);
  GenerateEdges();

  cout << "Graph Reading Done..." << endl;
}

// Graph &Graph::operator=(const Graph &other) {
//   size = other.size;
//   adjacency_table = other.adjacency_table;
//   vertex_to_aggregate = other.vertex_to_aggregate;
//   edges = other.edges;
//   edge_to_index = other.edge_to_index;
//   aggregates = other.aggregates;
//   interfaces = other.interfaces;
//   return *this;
// }

void Graph::GenerateEdges() {
  for (VertexIndex i = 0; i < adjacency_table.size(); ++i)
    for (auto j : adjacency_table[i]) {
      if (j > i) {
        Edge e(i, j);
        edges.push_back(e);
        edge_to_index[e] = edges.size() - 1;
      }
    }
}

void Graph::GenerateVertexMap() {
  for (AggregateIndex i = 0; i < aggregates.size(); ++i)
    for (auto vi : aggregates[i]->vertices)
      vertex_to_aggregate[vi] = i;
}

void Graph::UpdateVertexMap(const vector<AggregateIndex> &agg_indices) {
  for (auto i : agg_indices)
    for (auto vi : aggregates[i]->vertices)
      vertex_to_aggregate[vi] = i;
}

void Graph::GenerateInterfaces(AggregateIndex i) {
  for (auto edge : aggregates[i]->interface_edges) {
    AggregateIndex j = vertex_to_aggregate[edge.tail];
    interfaces.Add(i, j, edge);
  }
}

// void Graph::GeneratePiHColumns(const AggregatePtr &agg) {
//   for (auto edge : agg->interface_edges) {
//     if (agg->PiHColumns.count(edge.head)) continue;
//     unordered_set<VertexIndex> visited;
//     visited.insert(edge.head);
//     TreeSize(edge.head, agg, visited, agg->PiHColumns[edge.head]);
//   }
// }
//
// int Graph::TreeSize(VertexIndex root, const AggregatePtr &agg,
//                     unordered_set<VertexIndex> &visited, EdgeVector &v) {
//   vector<VertexIndex> children;
//   for (auto nbr : adjacency_table[root])
//     if (agg->HasVertex(nbr) && !visited.count(nbr)) {
//       children.push_back(nbr);
//       visited.insert(nbr);
//     }
//   int size = 0;
//   for (auto child : children) {
//     int subtreesize = TreeSize(child, agg, visited, v);
//     Edge e(root, child);
//     v.edge_indices.push_back(edge_to_index[e]);
//     double sign = child < root ? 1.0 : -1.0;
//     v.values.push_back(sign * (double)subtreesize / agg->Size());
//     size += subtreesize;
//   }
//   return size + 1;
// }

void Graph::DoMatching() {
  if (reshaped) {
    throw runtime_error(
        "The graph has been reshaped, no matching is performed!");
    return;
  }

  size_t n = NumOfAggregates();
  if (n == 1) {
    throw runtime_error("Only 1 aggregate, no matching is performed!");
    return;
  }

  // Table that stores the number of connections between aggregates
  struct ConnectivityEntry {
    AggregateIndex nbr_index;
    int num_connections;

    ConnectivityEntry(AggregateIndex i, int n):
        nbr_index(i), num_connections(n) {}
  };
  struct GoBefore {
    bool operator() (const ConnectivityEntry& e1, const ConnectivityEntry& e2)
        const {
      return e1.nbr_index < e2.nbr_index;
    }
  };
  std::vector<std::vector<ConnectivityEntry> > connectivity_table(n);

  for (auto interface : interfaces) {
    AggregateIndex i = interface.head_index, j = interface.tail_index;
    int size = interface.Size();
    connectivity_table[i].push_back(ConnectivityEntry(j, size));;
    connectivity_table[j].push_back(ConnectivityEntry(i, size));
  }
  for (auto &row : connectivity_table)
    sort(row.begin(), row.end(), GoBefore());

  vector<int> grouping_label(n, -1);
  vector<pair<int, AggregateIndex> > degree_index(n);
  for (AggregateIndex i = 0; i < n; ++i) {
    degree_index[i] = make_pair(aggregates[i]->interface_edges.size(), i);
  }
  sort(degree_index.begin(), degree_index.end());
  vector<AggregatePtr> aggregates_new;
  size_t count = 0;
  for (auto p : degree_index) {
    AggregateIndex i = p.second;
    if (grouping_label[i] != -1) continue;
    int max = 0, j = -1;
    for (auto entry : connectivity_table[i]) {
      if (grouping_label[entry.nbr_index] == -1 &&
          entry.num_connections > max) {
        max = entry.num_connections;
        j = entry.nbr_index;
      }
    }
    // If all neighbors are already matched
    if (j == -1) {
      unordered_map<AggregateIndex, int> num_connections_to;
      for (auto entry : connectivity_table[i]) {
        num_connections_to[grouping_label[entry.nbr_index]]
            += entry.num_connections;
      }
      ConnectivityEntry chosen(count, 0);
      for (auto entry : num_connections_to) {
        if (entry.second > chosen.num_connections ||
            (entry.second == chosen.num_connections &&
             entry.first < chosen.nbr_index)) {
          chosen = ConnectivityEntry(entry.first, entry.second);
        }
      }
      if (num_connections_to.empty()) {
        cout << "Graph not connected. Aborting!" << endl;
        exit(1);
      }
      grouping_label[i] = chosen.nbr_index;
      aggregates_new[chosen.nbr_index]->Merge(aggregates[i]);
    }
    else {
      aggregates_new.push_back(Merge(aggregates[i], aggregates[j]));
      grouping_label[i] = grouping_label[j] = count++;
    }
  }
  aggregates = aggregates_new;

  // Recompute the set of interfaces
  InterfaceSet interfaces_new;
  for (auto interface : interfaces) {
    AggregateIndex ai = grouping_label[interface.head_index],
                   aj = grouping_label[interface.tail_index];
    if (ai != aj) interfaces_new.Merge(ai, aj, interface.edges);
  }
  interfaces = interfaces_new;

  GenerateVertexMap();
}

void Graph::DoMatching(Graph *c_graph) {
  DoMatching();

  auto size = aggregates.size();
  for (VertexIndex i = 0; i < size; ++i) {
    AggregatePtr agg = make_shared<Aggregate>();
    agg->AddVertex(i);
    c_graph->aggregates.push_back(agg);
  }

  c_graph->adjacency_table.resize(size);
  for (auto interface : interfaces) {
    AggregateIndex i = interface.head_index, j = interface.tail_index;
    c_graph->adjacency_table[i].push_back(VertexIndex(j));
    c_graph->aggregates[i]->AddInterfaceEdge(Edge(i, j));
    c_graph->adjacency_table[j].push_back(VertexIndex(i));
    c_graph->aggregates[j]->AddInterfaceEdge(Edge(j, i));
    c_graph->interfaces.Add(i, j, Edge(i, j));
  }

  for (VertexIndex i = 0; i < size; ++i) {
    sort(c_graph->adjacency_table[i].begin(), c_graph->adjacency_table[i].end());
  }

  c_graph->size = size;
  c_graph->vertex_to_aggregate.resize(size);
  c_graph->GenerateEdges();
}

// void Graph::SplitAggregates(const std::vector<AggregateIndex> &indices) {
//   int size = indices.size();
//   for (auto i : indices) {
//     AggregatePtr a = aggregates[i];
//     size_t size;
//     while ((size = a->sub_aggregates.size()) == 1)
//       a = a->sub_aggregates[0];
//     if (size == 0) {
//       cout << "Aggregate " << i << " is a vertex! Skipped" << endl;
//       continue;
//     }
//
//     vector<AggregateIndex> agg_indices(size);
//     agg_indices[0] = i;
//     for (size_t j = 1, offset = aggregates.size() - 1; j < size; ++j)
//       agg_indices.push_back(offset + j);
//
//     interfaces.remove(i);
//     aggregates[i] = a->sub_aggregates[0];
//     aggregates.insert(aggregates.end(), a->sub_aggregates.begin()+1,
//                       a->sub_aggregates.end());
//
//     UpdateVertexMap(agg_indices);
//   }
//
//   // Generate new interfaces
//   for (auto i : indices) GenerateInterfaces(i);
//   for (size_t i = size; i < aggregates.size(); ++i) GenerateInterfaces(i);
// }

void Graph::GetAggregate(AggregateIndex i, std::vector<int> *vertices) const {
  for (auto v : aggregates[i]->vertices) {
    vertices->push_back(v);
  }
}

dCSRmat *Graph::GetWeightedLaplacian() const {
  size_t n = NumOfAggregates();

  INT nnz = n;
  for (AggregateIndex i = 0; i < n; ++i) {
    nnz += aggregates[i]->interface_edges.size();
  }

  dCOOmat *B = (dCOOmat*)malloc(sizeof(dCOOmat));
  B->row = n;
  B->col = n;
  B->nnz = nnz;
  B->rowind = (INT*)malloc(sizeof(INT)*nnz);
  B->colind = (INT*)malloc(sizeof(INT)*nnz);
  B->val = (REAL*)malloc(sizeof(REAL)*nnz);

  int ind = 0;
  for (AggregateIndex i = 0; i < n; ++i) {
    B->rowind[ind] = i;
    B->colind[ind] = i;
    B->val[ind] = (REAL)(aggregates[i]->interface_edges.size());
    ++ind;
  }
  for (auto interface : interfaces) {
    AggregateIndex i = interface.head_index, j = interface.tail_index;
    B->rowind[ind] = i;
    B->colind[ind] = j;
    B->val[ind] = -(REAL)interface.Size();
    ++ind;
    B->rowind[ind] = j;
    B->colind[ind] = i;
    B->val[ind] = -(REAL)interface.Size();
    ++ind;
  }

  dCSRmat *A = (dCSRmat*)malloc(sizeof(dCSRmat));
  dcoo_2_dcsr(B, A);
  free(B);

  return A;
}

// mfem::Vector Graph::GetCoarseRHS(const mfem::Vector &f) const {
//   assert(f.Size() == size);
//   size_t n = NumOfAggregates();
//   mfem::Vector fh(n);
//   for (AggregateIndex i = 0; i < n; ++i) {
//     fh(i) = 0;
//     for (VertexIndex vi : aggregates[i]->vertices)
//       fh(i) += f(vi);
//   }
//   return fh;
// }
//
// mfem::Vector Graph::Prolongate(const mfem::Vector &uh_c) const {
//   size_t n = Size(), nc = NumOfAggregates();
//   mfem::Vector uh(n);
//   for (AggregateIndex i = 0; i < nc; ++i)
//     for (VertexIndex vi : aggregates[i]->vertices)
//       uh(vi) = uh_c(i);
//   return uh;
// }
//
// mfem::SparseMatrix *Graph::GetCoarseEdgeSpaceBasis() {
//   assert(edge_to_index.size() == edges.size());
//   size_t ne = edge_to_index.size(), nec = interfaces.size();
//   mfem::SparseMatrix Phi_t(nec, ne);
//   int row = 0;
//   for (auto interface : interfaces) {
//     AggregatePtr head_agg = aggregates[interface.head_index],
//                  tail_agg = aggregates[interface.tail_index];
//     if (!head_agg->PiHGenerated()) GeneratePiHColumns(head_agg);
//     if (!tail_agg->PiHGenerated()) GeneratePiHColumns(tail_agg);
//     for (auto edge : interface.edges) {
//       VertexIndex head = edge.head, tail = edge.tail;
//       if (!head_agg->HasVertex(head)) swap(head, tail);
//       mfem::Array<int> cols1;
//       toArray(head_agg->PiHColumns[head].edge_indices, cols1);
//       Phi_t.AddRow(row, cols1, toVector(head_agg->PiHColumns[head].values));
//       mfem::Array<int> cols2;
//       toArray(tail_agg->PiHColumns[tail].edge_indices, cols2);
//       mfem::Vector values = toVector(tail_agg->PiHColumns[tail].values);
//       Phi_t.AddRow(-row-1, cols2, values);
//
//       Phi_t.Set(row, edge_to_index[edge], head < tail ? 1.0 : -1.0);
//     }
//     ++row;
//   }
//   Phi_t.Finalize();
//   Phi_t.SortColumnIndices();
//   return Transpose(Phi_t);
// }
//
// mfem::Vector Graph::GetLocalEstimates(mfem::Vector f, const mfem::Vector &uh,
//                                       mfem::Vector y, double CP) const {
//   G->AddMult(uh, y, -1.0);
//   G->AddMultTranspose(y, f, -1.0);
//   size_t nc = NumOfAggregates();
//   mfem::Vector local_est(nc);
//   local_est = 0.0;
//   // Vertices
//   for (AggregateIndex i = 0; i < nc; ++i)
//     for (VertexIndex vi : aggregates[i]->vertices)
//       local_est(i) += f(vi) * f(vi) / CP / CP;
//   // Interface edges
//   for (auto interface : interfaces) {
//     AggregateIndex i = interface.head_index, j = interface.tail_index;
//     for (auto edge : interface.edges) {
//       auto it = edge_to_index.find(edge);
//       assert(it != edge_to_index.end());
//       EdgeIndex ei = it->second;
//       double v = 0.5 * y(ei) * y(ei);
//       local_est(i) += v;
//       local_est(j) += v;
//     }
//   }
//   // Interior edges
//   for (auto pair : edge_to_index) {
//     const Edge &edge = pair.first;
//     EdgeIndex ei = pair.second;
//     VertexIndex vi = edge.head, vj = edge.tail;
//     AggregateIndex i;
//     if ((i = vertex_to_aggregate[vi]) == vertex_to_aggregate[vj])
//        local_est(i) += y(ei) * y(ei);
//   }
//   return local_est;
// }

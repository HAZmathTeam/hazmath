#include <algorithm>
#include <cassert>
#include <climits>
#include <fstream>
#include <iostream>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <unordered_map>

#include "algorithm.h"
#include "graph.h"

using namespace std;

/* Assumption on the mtx file representing a matrix:
 * 1. Vertex numbers start from 1 instead of 0
 * 2. For each edge (i,j) (i < j), the entry (j,i) is not recorded
 * 3. Neighbors of a vertex appear in ascending order
 */
Graph::Graph(const char *filename) {
  ifstream file(filename);
  string line;
  while (getline(file, line) && (line.empty() || line[0] == '%'))
    ;
  int nrows, ncols, nnz;
  istringstream iss(line);
  if (!(iss >> nrows >> ncols >> nnz)) {
    throw runtime_error("Error reading matrix sizes! Exiting...");
  }
  if (nrows != ncols) {
    throw logic_error("Not square matrix! Exiting...");
  }
  vector<vector<int>> adjacency_table(nrows);
  int count = 0, edge_count = 0, i, j;
  while (getline(file, line)) {
    ++count;
    istringstream iss(line);
    if (!(iss >> i >> j)) {
      throw runtime_error("Error reading matrix entry! Exiting...");
    }
    if (i == j)
      continue;
    --i;
    --j;
    assert(adjacency_table[i].empty() || adjacency_table[i].back() < j);
    adjacency_table[i].push_back(j);
    assert(adjacency_table[j].empty() || adjacency_table[j].back() < i);
    adjacency_table[j].push_back(i);
    edge_count += 2;
  }
  assert(count == nnz);
  file.close();

  A = (iCSRmat *)malloc(sizeof(iCSRmat));
  *A = icsr_create(nrows, ncols, edge_count);
  int ind = 0;
  for (int i = 0; i < nrows; ++i) {
    A->IA[i] = ind;
    for (auto j : adjacency_table[i]) {
      A->JA[ind] = j;
      A->val[ind] = 1;
      ++ind;
    }
  }
  A->IA[nrows] = ind;
  assert(ind == A->nnz);

  cout << "Graph Reading Done..." << endl;
}

std::vector<int> Graph::getNeighbors(int i) const {
  return std::vector<int>(A->JA + A->IA[i], A->JA + A->IA[i + 1]);
}

Graph::Graph(const Graph &other) {
  A = (iCSRmat *)malloc(sizeof(iCSRmat));
  *A = icsr_create(other.A->row, other.A->col, other.A->nnz);
  icsr_cp(other.A, A);
  aggregates = other.aggregates;
}

Graph &Graph::operator=(Graph other) {
  std::swap(A, other.A);
  std::swap(aggregates, other.aggregates);
  return *this;
}

void Graph::doMatching(Graph *c_graph) {
  int n = size();
  if (n == 1) {
    throw runtime_error("Only 1 node, no matching is performed!");
    return;
  }

  // Table that stores the number of connections between aggregates
  struct ConnectivityEntry {
    int nbr_index;
    int num_connections;

    ConnectivityEntry(int i, int n) : nbr_index(i), num_connections(n) {}
  };
  struct GoBefore {
    bool operator()(const ConnectivityEntry &e1,
                    const ConnectivityEntry &e2) const {
      return e1.nbr_index < e2.nbr_index;
    }
  };

  vector<vector<ConnectivityEntry>> connectivity_table(n);
  for (int i = 0; i < A->row; ++i) {
    for (int ind = A->IA[i]; ind < A->IA[i + 1]; ++ind) {
      int j = A->JA[ind];
      int w = A->val[ind];
      connectivity_table[i].push_back(ConnectivityEntry(j, w));
      connectivity_table[j].push_back(ConnectivityEntry(i, w));
    }
  }
  for (auto &row : connectivity_table)
    sort(row.begin(), row.end(), GoBefore());

  vector<int> grouping_label(n, -1);
  vector<pair<int, int>> degree_index(n);
  for (int i = 0; i < n; ++i) {
    int degree = 0;
    for (auto entry : connectivity_table[i]) {
      degree += entry.num_connections;
    }
    degree_index[i] = make_pair(degree, i);
  }
  sort(degree_index.begin(), degree_index.end());
  int count = 0;
  for (auto p : degree_index) {
    int i = p.second;
    if (grouping_label[i] != -1)
      continue;
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
      unordered_map<int, int> num_connections_to;
      for (auto entry : connectivity_table[i]) {
        num_connections_to[grouping_label[entry.nbr_index]] +=
            entry.num_connections;
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
      aggregates[chosen.nbr_index].push_back(i);
    } else {
      aggregates.push_back(vector<int>{i, j});
      grouping_label[i] = grouping_label[j] = count++;
    }
  }

  vector<set<int>> c_adjacency_table(count);
  for (int i = 0; i < A->row; ++i) {
    for (int ind = A->IA[i]; ind < A->IA[i + 1]; ++ind) {
      int j = A->JA[ind];
      int I = grouping_label[i];
      int J = grouping_label[j];
      if (I != J) {
        c_adjacency_table[I].insert(J);
      }
    }
  }
  int nnz = 0;
  for (auto nbs : c_adjacency_table) {
    nnz += nbs.size();
  }
  iCSRmat *c_A = (iCSRmat *)malloc(sizeof(iCSRmat));
  *c_A = icsr_create(count, count, nnz);

  int ind = 0;
  for (int i = 0; i < count; ++i) {
    c_A->IA[i] = ind;
    for (auto nb : c_adjacency_table[i]) {
      c_A->JA[ind] = nb;
      c_A->val[ind] = 1;
      ++ind;
    }
  }
  c_A->IA[count] = ind;
  assert(ind == nnz);

  c_graph->A = c_A;
}

void Graph::doMatchingDegreeBased(Graph *c_graph, int seed) {
  int n = size();
  if (n == 1) {
    throw runtime_error("Only 1 node, no matching is performed!");
    return;
  }

  vector<int> degrees(n, 0);
  for (int i = 0; i < A->row; ++i) {
    for (int ind = A->IA[i]; ind < A->IA[i + 1]; ++ind) {
      int j = A->JA[ind];
      int w = A->val[ind];
      degrees[i] += w;
      degrees[j] += w;
    }
  }

  vector<int> vertices(n);
  std::iota(vertices.begin(), vertices.end(), 0);
  // Randomly shuffle vertices.
  std::shuffle(vertices.begin(), vertices.end(),
               std::default_random_engine(seed));

  // Set up uniform distribution on [0, 1).
  std::mt19937 gen(seed);
  std::uniform_real_distribution<> dis(0.0, 1.0);

  auto pickMinUniformly = [&](int value, int &min, int &num_min) -> bool {
    if (value == min) {
      if (dis(gen) < 1.0 / ++num_min) {
        return true;
      }
    } else if (value < min) {
      num_min = 1;
      min = value;
      return true;
    }
    return false;
  };

  vector<int> grouping_label(n, -1);
  vector<int> cardinalities;
  int agg_count = 0;
  for (auto i : vertices) {
    if (grouping_label[i] != -1) {
      continue;
    }

    // Randomly pick the neighboring vertex with smallest degree.
    int selected_nbr = -1;
    int min = INT_MAX;
    int num_min = 0;
    std::set<int> nbr_aggregates;
    for (int ind = A->IA[i]; ind < A->IA[i + 1]; ++ind) {
      int j = A->JA[ind];
      if (grouping_label[j] == -1) {
        if (pickMinUniformly(degrees[j], min, num_min)) {
          selected_nbr = j;
        }
      } else {
        nbr_aggregates.insert(grouping_label[j]);
      }
    }

    // If all neighbors are already matched.
    if (selected_nbr == -1) {
      // Randomly pick the neighboring aggregate with smallest cardinality.
      int selected_agg;
      int min = INT_MAX;
      int num_min = 0;
      for (auto agg : nbr_aggregates) {
        if (pickMinUniformly(cardinalities[agg], min, num_min)) {
          selected_agg = agg;
        }
      }
      grouping_label[i] = selected_agg;
      aggregates[selected_agg].push_back(i);
      ++cardinalities[selected_agg];
    } else {
      aggregates.push_back(vector<int>{i, selected_nbr});
      grouping_label[i] = agg_count;
      grouping_label[selected_nbr] = agg_count;
      --degrees[i];
      --degrees[selected_nbr];
      cardinalities.push_back(2);
      ++agg_count;
    }
  }

  vector<set<int>> c_adjacency_table(agg_count);
  for (int i = 0; i < A->row; ++i) {
    for (int ind = A->IA[i]; ind < A->IA[i + 1]; ++ind) {
      int j = A->JA[ind];
      int I = grouping_label[i];
      int J = grouping_label[j];
      if (I != J) {
        c_adjacency_table[I].insert(J);
      }
    }
  }
  int nnz = 0;
  for (auto nbs : c_adjacency_table) {
    nnz += nbs.size();
  }
  iCSRmat *c_A = (iCSRmat *)malloc(sizeof(iCSRmat));
  *c_A = icsr_create(agg_count, agg_count, nnz);

  int ind = 0;
  for (int i = 0; i < agg_count; ++i) {
    c_A->IA[i] = ind;
    for (auto nb : c_adjacency_table[i]) {
      c_A->JA[ind] = nb;
      c_A->val[ind] = 1;
      ++ind;
    }
  }
  c_A->IA[agg_count] = ind;
  assert(ind == nnz);

  c_graph->A = c_A;
}

void Graph::getAggregate(int i, std::vector<int> *vertices) const {
  for (auto v : aggregates[i]) {
    vertices->push_back(v);
  }
}

dCSRmat *Graph::getWeightedLaplacian() const {
  int n = size();
  int nnz = A->nnz + n;

  dCSRmat *L = (dCSRmat *)malloc(sizeof(dCSRmat));
  *L = dcsr_create(n, n, nnz);

  int L_ind = 0;
  for (int i = 0; i < n; ++i) {
    L->IA[i] = L_ind;
    REAL sum = 0;
    int L_ind_diag;
    int ind = A->IA[i];
    while (ind < A->IA[i + 1]) {
      int j = A->JA[ind];
      assert(j != i);
      if (j > i) {
        break;
      }
      L->JA[L_ind] = j;
      L->val[L_ind] = -A->val[ind];
      sum += A->val[ind];
      ++L_ind;
      ++ind;
    }
    L_ind_diag = L_ind++;
    while (ind < A->IA[i + 1]) {
      int j = A->JA[ind];
      assert(j != i);
      L->JA[L_ind] = j;
      L->val[L_ind] = -A->val[ind];
      sum += A->val[ind];
      ++L_ind;
      ++ind;
    }
    L->JA[L_ind_diag] = i;
    L->val[L_ind_diag] = sum;
  }
  L->IA[n] = L_ind;
  assert(L_ind == nnz);
  return L;
}

vector<int> Graph::getHamiltonianPath(int seed) const {
  srand(seed);
  INT root = rand() % size();

  // Build spanning tree.
  Tree *tree = new Tree(root);
  queue<Tree *> q;
  q.push(tree);
  vector<bool> discovered(size(), false);
  discovered[root] = true;
  while (!q.empty()) {
    auto curr = q.front();
    for (auto neighbor : getNeighbors(curr->vertex)) {
      if (!discovered[neighbor]) {
        curr->children.push_back(new Tree(neighbor));
        q.push(curr->children.back());
        discovered[neighbor] = true;
      }
    }
    q.pop();
  }

  return ::getHamiltonianPath(tree);
}

#include "algorithm.h"
#include "graph.h"
#include "hamiltonian_algorithm.h"
#include <algorithm>
#include <functional>
#include <iostream>
#include <numeric>



extern "C" {
void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda, double *w,
            double *work, int *lwork, int *info);
}

// std::string genFilePrefix(const std::string &filename,
//                           const std::string &matching_alg,
//                           const std::string &compression_alg) {
//   return filename.substr(0, filename.find_last_of(".")) + "_" + matching_alg
//   +
//          "_" + compression_alg + "_";
// }

void AggregationBasedAlgorithm::setupHierarchy(Graph graph,
                                               std::vector<dCSRmat *> &Qj_array,
                                               std::vector<int> &Nj_array,
                                               bool weighted) const {
  int n = graph.size();

  /*
  dCSRmat *Q = (dCSRmat*)malloc(sizeof(dCSRmat));
  dCSRmat identity = dcsr_create_identity_matrix(n, 0);
  dcsr_cp(&identity, Q);
  */

  int level = 0;
  Nj_array.push_back(n);
  // Recursively perform aggregation algorithm until there is only one vertex in
  // the graph.
  while (graph.size() > 1) {
    if (print_level == 0) {
      std::cout << "Graph size: " << graph.size() << std::endl;
    }
    dCSRmat *A = graph.getLaplacian();

    Graph c_graph;
    doMatching(graph, &c_graph);
    assert(graph.numOfAggregates() == c_graph.size());
    int Nj = graph.size();
    int nj = c_graph.size();
    Nj_array.push_back(nj);
    int numBlocks = getNumBlocks(1 << level);

    int Qj_nnz = 0;
    for (int i = 0; i < nj; ++i) {
      int size = graph.getAggregateSize(i);
      Qj_nnz += numBlocks * size * size;
    }
    Qj_nnz += n - Nj * numBlocks;
    dCOOmat *Qj_coo = (dCOOmat *)malloc(sizeof(dCOOmat));
    Qj_coo->row = n;
    Qj_coo->col = n;
    Qj_coo->nnz = Qj_nnz;
    Qj_coo->rowind = (INT *)malloc(sizeof(INT) * Qj_nnz);
    Qj_coo->colind = (INT *)malloc(sizeof(INT) * Qj_nnz);
    Qj_coo->val = (REAL *)malloc(sizeof(REAL) * Qj_nnz);
    INT Qj_coo_ind = 0;

    // Iterate over the aggregates in the coarse graph.
    for (int i = 0, count = 0; i < nj; ++i) {
      std::vector<int> vertices;
      graph.getAggregate(i, &vertices);
      int ni = vertices.size();
      sort(vertices.begin(), vertices.end());

      dCSRmat *Ai = (dCSRmat *)malloc(sizeof(dCSRmat));
      dcsr_getblk(A, vertices.data(), vertices.data(), ni, ni, Ai);

      double a[ni][ni];
      memset(*a, 0, ni * ni * sizeof(REAL));
      for (INT i = 0; i < Ai->row; ++i) {
        double rowsum = 0.0;
        for (INT ind = Ai->IA[i]; ind < Ai->IA[i + 1]; ++ind) {
          INT j = Ai->JA[ind];
          if (j != i) {
            a[i][j] = Ai->val[ind];
            if (weighted) {
              a[i][j] /= graph.getAbsCard(i) * graph.getAbsCard(j);
            };
            rowsum += Ai->val[ind];
          }
        }
        a[i][i] = -rowsum;
        if (weighted) {
          a[i][i] /= graph.getAbsCard(i) * graph.getAbsCard(i);
        };
      }
      // dcsr_write_dcoo("Ai.dat", Ai);
      // print_full_mat(ni, ni, *a, "a");

      char jobz = 'V', uplo = 'U';
      double w[ni];
      int n5 = 5 * ni;
      double *work = (double *)calloc(n5, sizeof(double));
      int info;
      dsyev_(&jobz, &uplo, &ni, *a, &ni, w, work, &n5, &info);
      free(work);
      if (info) {
        std::cout << "Eigenvalue computations error; Error code: " << info << std::endl;
        return;
      }

      for (int l = 0; l < numBlocks; ++l) {
        int k = 0;
        while (k < 2) {
          for (int ind = 0; ind < ni; ++ind) {
            Qj_coo->rowind[Qj_coo_ind] = 2 * nj * l + k * nj + i;
            Qj_coo->colind[Qj_coo_ind] = vertices[ind];
            Qj_coo->val[Qj_coo_ind] = a[ind][k];
            // if (isnan(a[ind][k]))
            //   std::cout << "Found NaN!" << std::endl;
            ++Qj_coo_ind;
          }
          // std::cout << "Row: " << 2 * nj * l + k * nj + i << std::endl;
          // std::cout << "Values: " << v(0) << " " << v(1) << " " << v.Norml2()
          //      << std::endl;
          ++k;
        }
        while (k < ni) {
          for (int ind = 0; ind < ni; ++ind) {
            Qj_coo->rowind[Qj_coo_ind] =
                2 * nj * numBlocks + l * (Nj - 2 * nj) + count + (k - 2);
            Qj_coo->colind[Qj_coo_ind] = vertices[ind];
            Qj_coo->val[Qj_coo_ind] = a[ind][k];
            /* if (isnan(a[ind][k])) {
              std::cout << "Found NaN!" << std::endl
                   << "ni: " << ni << std::endl
                   << "k: " << k << std::endl;
            } */
            ++Qj_coo_ind;
          }
          // std::cout << "Row: "
          //      << 2 * nj * numBlocks + l * (Nj - 2 * nj) + count + (k - 2)
          //      << std::endl;
          ++k;
        }
        for (auto i = 0; i < ni; ++i) {
          vertices[i] += Nj;
        }

        // std::cout << l << std::endl;
      }
      count += ni - 2;
    }

    for (int i = Nj * numBlocks; i < n; ++i) {
      Qj_coo->rowind[Qj_coo_ind] = i;
      Qj_coo->colind[Qj_coo_ind] = i;
      Qj_coo->val[Qj_coo_ind] = 1.0;
      ++Qj_coo_ind;
    }

    dCSRmat *Qj = (dCSRmat *)malloc(sizeof(dCSRmat));
    dcoo_2_dcsr(Qj_coo, Qj);
    free(Qj_coo->rowind);
    free(Qj_coo->colind);
    free(Qj_coo->val);
    free(Qj_coo);

    /*
    dCSRmat *Q1 = (dCSRmat*)malloc(sizeof(dCSRmat));
    dcsr_mxm(Qj, Q, Q1);
    dcsr_free(Q);
    // dcsr_free(Qj);
    Q = Q1;
    */

    Qj_array.push_back(Qj);

    graph = std::move(c_graph);
    ++level;
    dcsr_free(A);
  }
}

// void AggregationBasedAlgorithm::readStoredHierarchy(
//     Graph graph, std::string matching_alg, std::string compression_alg,
//     std::vector<dCSRmat *> &Qj_array, std::vector<int> &Nj_array) const {
//   const std::string prefix =
//       genFilePrefix(graph.getFilename(), matching_alg, compression_alg);
//   int n = graph.size();
// }

// void AggregationBasedAlgorithm::saveHierarchy(
//     Graph graph, std::string matching_alg, std::string compression_alg,
//     std::vector<dCSRmat *> &Qj_array, std::vector<int> &Nj_array) const {
//   const std::string prefix =
//       genFilePrefix(graph.getFilename(), matching_alg, compression_alg);
//   const int n = graph.size();
// }

REAL AggregationBasedAlgorithm::compAndDecomp(int n, double *v,
                                              const std::vector<dCSRmat *> &Qj_array,
                                              int largestK, double *v2) const {
  /*
  REAL vt[n];
  dcsr_mxv(Q, v, vt);

  std::vector<double> vt_sort(n);
  for (int i = 0; i < n; ++i) {
    vt_sort[i] = abs(vt[i]);
  }
  sort(vt_sort.begin(), vt_sort.end(), std::greater<double>());
  for (int i = 0; i < n; ++i) {
    if (abs(vt[i]) < vt_sort[largestK]) {
      vt[i] = 0;
    }
  }

  dCSRmat *Qt = (dCSRmat*)malloc(sizeof(dCSRmat));
  dcsr_trans(Q, Qt);
  // dcsr_write_dcoo("Q.dat", Q);
  // dcsr_write_dcoo("Qt.dat", Qt);
  REAL v1[n];
  dcsr_mxv(Qt, vt, v1);
  REAL e[n];
  array_axpyz(n, -1.0, v, v1, e);

  for (int i = 0; i < Q->nnz; ++i) {
    if (isnan(Q->val[i])) std::cout << i << " " << Q->JA[i] << std::endl;
  }
  for (int i = 0; i < Qt->nnz; ++i) {
    if (isnan(Qt->val[i])) std::cout << i << " " << Qt->JA[i] << std::endl;
  }

  std::cout << "Norm of vector ||v||: " << array_norm2(n, v) << std::endl
       << "Norm of error  ||v-v1||: " << array_norm2(n, e) << std::endl
       << "Relative error ||v-v1||/||v||:  "
       << array_norm2(n, e) / array_norm2(n, v) << std::endl;

  dcsr_free(Q);
  dcsr_free(Qt);
  */

  REAL vj[n];
  array_cp(n, v, vj);
  for (auto Qj : Qj_array) {
    REAL v_temp[n];
    dcsr_mxv(Qj, vj, v_temp);
    array_cp(n, v_temp, vj);
  }

  std::vector<double> vt_sort(n);
  for (int i = 0; i < n; ++i) {
    vt_sort[i] = abs(vj[i]);
  }
  sort(vt_sort.begin(), vt_sort.end(), std::greater<double>());
  for (int i = 0; i < n; ++i) {
    if (abs(vj[i]) < vt_sort[largestK]) {
      vj[i] = 0;
    }
  }

  array_cp(n, vj, v2);
  for (auto it = Qj_array.rbegin(); it != Qj_array.rend(); ++it) {
    dCSRmat *Qj_t = (dCSRmat *)malloc(sizeof(dCSRmat));
    dcsr_trans(*it, Qj_t);
    REAL v_temp[n];
    dcsr_mxv(Qj_t, v2, v_temp);
    array_cp(n, v_temp, v2);
    dcsr_free(Qj_t);
  }
  REAL e2[n];
  array_axpyz(n, -1.0, v, v2, e2);

  if (print_level == 0) {
    std::cout << std::endl
         << "Plain Encoding" << std::endl
         << "Norm of vector ||v||: " << array_norm2(n, v) << std::endl
         << "Norm of error  ||v-v2||: " << array_norm2(n, e2) << std::endl
         << "Relative error ||v-v2||/||v||: "
         << array_norm2(n, e2) / array_norm2(n, v) << std::endl;
  }

  return array_norm2(n, e2) / array_norm2(n, v);
}

void AggregationBasedAlgorithm::compAndDecompAdaptive(
    int n, double *v, const std::vector<dCSRmat *> &Qj_array,
    const std::vector<int> &Nj_array, int largestK, double p, double *v3) const {
  /* ------------------ Testing adaptive encoding -------------------- */
  std::vector<std::vector<REAL>> vj_array{std::vector<REAL>(v, v + n)};
  for (auto Qj : Qj_array) {
    vj_array.push_back(std::vector<REAL>(n));
    dcsr_mxv(Qj, (vj_array.rbegin() + 1)->data(), vj_array.back().data());
  }
  // Finding the optimal basis
  assert(Nj_array.size() == vj_array.size());
  int num_levels = Nj_array.size();
  std::vector<std::vector<REAL>> sums(num_levels);
  std::vector<std::vector<bool>> labels(num_levels);
  for (int j = num_levels - 1; j >= 0; --j) {
    int i = 0, aux, offset;
    if (j < num_levels - 1) {
      aux = Nj_array[j] - 2 * Nj_array[j + 1];
      offset = (1 << (j + 1)) * Nj_array[j + 1];
    }
    for (int k = 0; k < (1 << j); ++k) {
      sums[j].push_back(0.0);
      for (int l = 0; l < Nj_array[j]; ++l) {
        sums[j].back() += pow(abs(vj_array[j][i++]), p);
      }
      if (j < num_levels - 1) {
        REAL sum_aux =
            accumulate(vj_array[j + 1].begin() + offset + k * aux,
                       vj_array[j + 1].begin() + offset + (k + 1) * aux, 0);
        REAL sum_children =
            sums[j + 1][2 * k] + sums[j + 1][2 * k + 1] + sum_aux;
        if (sum_children >= sums[j][k]) {
          labels[j].push_back(true);
        } else {
          labels[j].push_back(false);
          sums[j][k] = sum_children;
        }
      } else {
        labels[j].push_back(true);
      }
    }
  }
  // Encode
  std::vector<REAL> v_e;
  v_e.reserve(n);
  std::function<void(int, int)> encode = [&](int j, int k) {
    if (labels[j][k]) {
      v_e.insert(v_e.end(), vj_array[j].begin() + k * Nj_array[j],
                 vj_array[j].begin() + (k + 1) * Nj_array[j]);
    } else {
      encode(j + 1, 2 * k);
      encode(j + 1, 2 * k + 1);
      int offset = (1 << (j + 1)) * Nj_array[j + 1];
      int aux = Nj_array[j] - 2 * Nj_array[j + 1];
      v_e.insert(v_e.end(), vj_array[j + 1].begin() + offset + k * aux,
                 vj_array[j + 1].begin() + offset + (k + 1) * aux);
    }
  };
  encode(0, 0);
  // Trunk
  std::vector<double> vt_sort(n);
  for (int j = 0; j < n; ++j) {
    vt_sort[j] = abs(v_e[j]);
  }
  sort(vt_sort.begin(), vt_sort.end(), std::greater<double>());
  for (int j = 0; j < n; ++j) {
    if (abs(v_e[j]) < vt_sort[largestK]) {
      v_e[j] = 0;
    }
  }
  // Decompress
  auto it = v_e.begin();
  std::function<std::vector<REAL>(int, int)> decode = [&](int j, int k) -> std::vector<REAL> {
    if (labels[j][k]) {
      auto res = std::vector<REAL>(it, it + Nj_array[j]);
      it += Nj_array[j];
      return res;
    } else {
      auto v1 = decode(j + 1, 2 * k);
      auto v2 = decode(j + 1, 2 * k + 1);
      std::vector<REAL> v_segment;
      int Nj = Nj_array[j];
      v_segment.reserve(Nj);
      v_segment.insert(v_segment.end(), v1.begin(), v1.end());
      v_segment.insert(v_segment.end(), v2.begin(), v2.end());
      int aux = Nj - 2 * Nj_array[j + 1];
      v_segment.insert(v_segment.end(), it, it + aux);
      it += aux;
      std::vector<int> Is, Js;
      for (int i = 0; i < 2 * Nj_array[j + 1]; ++i) {
        Is.push_back(i);
        Js.push_back(i);
      }
      int offset = (1 << (j + 1)) * Nj_array[j + 1];
      for (int i = 0; i < aux; ++i) {
        Is.push_back(offset + i);
        Js.push_back(2 * Nj_array[j + 1] + i);
      }
      dCSRmat *Qj_block = (dCSRmat *)malloc(sizeof(dCSRmat));
      dcsr_getblk(Qj_array[j], Is.data(), Js.data(), Nj, Nj, Qj_block);
      dCSRmat *Qj_block_t = (dCSRmat *)malloc(sizeof(dCSRmat));
      dcsr_trans(Qj_block, Qj_block_t);
      std::vector<REAL> v_res(Nj);
      dcsr_mxv(Qj_block_t, v_segment.data(), v_res.data());
      return v_res;
    }
  };
  auto v3_vector = decode(0, 0);
  copy(v3_vector.begin(), v3_vector.end(), v3);
  // Compute the error
  REAL e3[n];
  array_axpyz(n, -1.0, v, v3, e3);
  if (print_level == 0) {
    std::cout << std::endl
         << "Adaptive Encoding" << std::endl
         << "Norm of vector ||v||: " << array_norm2(n, v) << std::endl
         << "Norm of error  ||v-v3||: " << array_norm2(n, e3) << std::endl
         << "Relative error ||v-v3||/||v||: "
         << array_norm2(n, e3) / array_norm2(n, v) << std::endl;
  }
}

void AggregationBasedAlgorithm::compAndDecomp(
    const Graph &graph, std::vector<REAL *> vectors,
    const std::vector<int> &nums_terms, const double p) const {
  std::cout << "Matching algorithm: " << matchingAlgorithm() << std::endl
            << "Compression algorithm: " << compressionAlgorithm() << std::endl;
  std::vector<dCSRmat *> Qj_array;
  std::vector<int> Nj_array;
  // try {
  //   // readStoredHierarchy(graph.getFilename(), matchingAlgorithm(),
  //   //                     compressionAlgorithm(), Qj_array, Nj_array);
  // } catch (exception &e) {
  // setupHierarchy(graph, Qj_array, Nj_array);
  // saveHierarchy(graph.getFilename(), matchingAlgorithm(),
  //               compressionAlgorithm(), Qj_array, Nj_array);
  // }
  setupHierarchy(graph, Qj_array, Nj_array);
  int n = graph.size();
  REAL *v2 = (REAL *)malloc(sizeof(REAL) * n);
  REAL *v3 = (REAL *)malloc(sizeof(REAL) * n);

  std::cout << "Average errors:" << std::endl;
  for (auto num_terms : nums_terms) {
    REAL error_sum = 0.0;
    for (REAL *v : vectors) {
      error_sum += compAndDecomp(n, v, Qj_array, num_terms, v2);
      if (isAdaptive()) {
        compAndDecompAdaptive(n, v, Qj_array, Nj_array, num_terms, p, v3);
      }
    }
    std::cout << num_terms << " terms: " << error_sum / vectors.size() << std::endl;
  }

  for (auto Qj : Qj_array) {
    dcsr_free(Qj);
  }
  free(v2);
  free(v3);
}

void HamiltonianAlgorithm::compAndDecomp(
    const Graph &graph, std::vector<REAL *> vectors,
    const std::vector<int> &nums_terms) const {
  std::cout << "Hamiltonian algorithm: " << std::endl;
  auto permutation = graph.getHamiltonianPath();
  int n = permutation.size();
  int N = 1;
  int L = 0;
  while (N < n) {
    N <<= 1;
    ++L;
  }

  const auto walsh_basis = getWalshBasis(L);
  const auto samples = getRandomProlongation(n, N, 0);
  std::cout << "Average errors:" << std::endl;
  for (auto num_terms : nums_terms) {
    REAL error_sum = 0.0;
    for (REAL *v : vectors) {
      const std::vector<REAL> prolongated_v =
          prolongate(std::vector<REAL>(v, v + n), permutation, samples);
      const std::vector<REAL> k_term_approximation =
          approximate(prolongated_v, walsh_basis, N, num_terms);
      std::vector<REAL> approximation =
          project(k_term_approximation, samples, permutation);

      REAL e2[n];
      array_axpyz(n, -1.0, v, approximation.data(), e2);

      if (print_level == 0) {
        std::cout << std::endl
             << "Norm of vector ||v||: " << array_norm2(n, v) << std::endl
             << "Norm of error  ||v-v2||: " << array_norm2(n, e2) << std::endl
             << "Relative error ||v-v2||/||v||: "
             << array_norm2(n, e2) / array_norm2(n, v) << std::endl;
      }

      error_sum += array_norm2(n, e2) / array_norm2(n, v);
    }
    std::cout << num_terms << " terms: " << error_sum / vectors.size() << std::endl;
  }

  deleteArray(walsh_basis);
}

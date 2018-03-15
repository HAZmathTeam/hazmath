#include <iostream>
#include <numeric>
#include <algorithm>
#include <functional>
#include "graph.hpp"

using namespace std;

extern "C" {
  void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda, double *w,
              double *work, int *lwork, int *info);
}

REAL *InitializeRHS(dCSRmat *A, int num_iterations = 100);

void setup_hierarchy(const char *file, dCSRmat *&A, vector<dCSRmat *> &Qj_array,
    vector<int> &Nj_array) {
  Graph graph(file);
  int n = graph.Size();
  A = graph.GetWeightedLaplacian();

  /*
  dCSRmat *Q = (dCSRmat*)malloc(sizeof(dCSRmat));
  dCSRmat identity = dcsr_create_identity_matrix(n, 0);
  dcsr_cp(&identity, Q);
  */

  int level = 0;
  Nj_array.push_back(n);
  while (graph.Size() > 1) {
    cout << "Graph size: " << graph.Size() << endl;
    dCSRmat *A = graph.GetWeightedLaplacian();

    Graph c_graph;
    graph.DoMatching(&c_graph);
    assert(graph.NumOfAggregates() == c_graph.Size());
    int Nj = graph.Size();
    int nj = c_graph.Size();
    Nj_array.push_back(nj);
    int numBlocks = 1 << level;

    int Qj_nnz = 0;
    for (int i = 0; i < nj; ++i) {
      int size = graph.GetAggregateSize(i);
      Qj_nnz += numBlocks * size * size;
    }
    Qj_nnz += n - Nj * numBlocks;
    dCOOmat *Qj_coo = (dCOOmat*)malloc(sizeof(dCOOmat));
    Qj_coo->row = n;
    Qj_coo->col = n;
    Qj_coo->nnz = Qj_nnz;
    Qj_coo->rowind = (INT*)malloc(sizeof(INT)*Qj_nnz);
    Qj_coo->colind = (INT*)malloc(sizeof(INT)*Qj_nnz);
    Qj_coo->val = (REAL*)malloc(sizeof(REAL)*Qj_nnz);
    INT Qj_coo_ind = 0;

    for (int i = 0, count = 0; i < nj; ++i) {
      vector<int> vertices;
      graph.GetAggregate(i, &vertices);
      int ni = vertices.size();
      sort(vertices.begin(), vertices.end());

      dCSRmat *Ai = (dCSRmat*)malloc(sizeof(dCSRmat));
      dcsr_getblk(A, vertices.data(), vertices.data(), ni, ni, Ai);

      double a[ni][ni];
      memset(*a, 0, ni*ni*sizeof(REAL));
      for (INT i = 0; i < Ai->row; ++i) {
        double rowsum = 0.0;
        for (INT ind = Ai->IA[i]; ind < Ai->IA[i+1]; ++ind) {
          INT j = Ai->JA[ind];
          if (j != i) {
            // cout << i << " "<< j << " " << Ai->val[ind] << endl;
            a[i][j] = Ai->val[ind];
            rowsum += Ai->val[ind];
          }
        }
        a[i][i] = -rowsum;
      }
      // dcsr_write_dcoo("Ai.dat", Ai);
      // print_full_mat(ni, ni, *a, "a");

      char jobz = 'V', uplo = 'U';
      double w[ni];
      int n5 = 5 * ni;
      double *work = (double *)calloc(n5,sizeof(double));
      int info;
      dsyev_(&jobz, &uplo, &ni, *a, &ni, w, work, &n5, &info);
      free(work);
      if (info) {
        cout << "Eigenvalue computations error; Error code: " << info << endl;
        return;
      }

      for (int l = 0; l < numBlocks; ++l) {
        int k = 0;
        while (k < 2) {
          for (int ind = 0; ind < ni; ++ind) {
            Qj_coo->rowind[Qj_coo_ind] = 2*nj*l + k*nj + i;
            Qj_coo->colind[Qj_coo_ind] = vertices[ind];
            Qj_coo->val[Qj_coo_ind] = a[ind][k];
            // if (isnan(a[ind][k])) cout << "Found NaN!" << endl;
            ++Qj_coo_ind;
          }
          // cout << "Row: " << 2*nj*l + k*nj + i << endl;
          // cout << "Values: " << v(0) << " " << v(1) << " " << v.Norml2() << endl;
          ++k;
        }
        while (k < ni) {
          for (int ind = 0; ind < ni; ++ind) {
            Qj_coo->rowind[Qj_coo_ind] = 2*nj*numBlocks + l*(Nj-2*nj) + count + (k-2);
            Qj_coo->colind[Qj_coo_ind] = vertices[ind];
            Qj_coo->val[Qj_coo_ind] = a[ind][k];
            /* if (isnan(a[ind][k])) {
              cout << "Found NaN!" << endl
                   << "ni: " << ni << endl
                   << "k: " << k << endl;
            } */
            ++Qj_coo_ind;
          }
          // cout << "Row: " << 2*nj*numBlocks + l*(Nj-2*nj) + count + (k-2) << endl;
          ++k;
        }
        for (auto i = 0; i < ni; ++i) {
          vertices[i] += Nj;
        }

        // cout << l << endl;
      }
      count += ni - 2;
    }

    for (int i = Nj*numBlocks; i < n; ++i) {
      Qj_coo->rowind[Qj_coo_ind] = i;
      Qj_coo->colind[Qj_coo_ind] = i;
      Qj_coo->val[Qj_coo_ind] = 1.0;
      ++Qj_coo_ind;
    }

    dCSRmat *Qj = (dCSRmat*)malloc(sizeof(dCSRmat));
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

    graph = c_graph;
    ++level;
    dcsr_free(A);
  }
}

void comp_decomp(double *v, dCSRmat *A, const vector<dCSRmat *> &Qj_array,
    const vector<int> &Nj_array, int threshold, double p, double* v2,
    double *v3) {
  int n = A->row;

  // Compress/decompress a smooth vector and compute the error
  if (v == NULL) {
    v = InitializeRHS(A);
  }

  /*
  REAL vt[n];
  dcsr_mxv(Q, v, vt);

  vector<double> vt_sort(n);
  for (int i = 0; i < n; ++i) {
    vt_sort[i] = abs(vt[i]);
  }
  sort(vt_sort.begin(), vt_sort.end(), greater<double>());
  for (int i = 0; i < n; ++i) {
    if (abs(vt[i]) < vt_sort[threshold]) {
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
    if (isnan(Q->val[i])) cout << i << " " << Q->JA[i] << endl;
  }
  for (int i = 0; i < Qt->nnz; ++i) {
    if (isnan(Qt->val[i])) cout << i << " " << Qt->JA[i] << endl;
  }

  cout << "Norm of vector ||v||: " << array_norm2(n, v) << endl
       << "Norm of error  ||v-v1||: " << array_norm2(n, e) << endl
       << "Relative error ||v-v1||/||v||:  "
       << array_norm2(n, e) / array_norm2(n, v) << endl;

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

  vector<double> vt_sort(n);
  for (int i = 0; i < n; ++i) {
    vt_sort[i] = abs(vj[i]);
  }
  sort(vt_sort.begin(), vt_sort.end(), greater<double>());
  for (int i = 0; i < n; ++i) {
    if (abs(vj[i]) < vt_sort[threshold]) {
      vj[i] = 0;
    }
  }

  array_cp(n, vj, v2);
  for (auto it = Qj_array.rbegin(); it != Qj_array.rend(); ++it) {
    dCSRmat *Qj_t = (dCSRmat*)malloc(sizeof(dCSRmat));
    dcsr_trans(*it, Qj_t);
    REAL v_temp[n];
    dcsr_mxv(Qj_t, v2, v_temp);
    array_cp(n, v_temp, v2);
    dcsr_free(Qj_t);
  }
  REAL e2[n];
  array_axpyz(n, -1.0, v, v2, e2);

  cout << endl << "Plain Encoding" << endl
       << "Norm of vector ||v||: " << array_norm2(n, v) << endl
       << "Norm of error  ||v-v2||: " << array_norm2(n, e2) << endl
       << "Relative error ||v-v2||/||v||: "
       << array_norm2(n, e2) / array_norm2(n, v) << endl;

  /* ------------------ Testing adaptive encoding -------------------- */
  vector<vector<REAL>> vj_array{vector<REAL>(v, v+n)};
  for (auto Qj : Qj_array) {
    vj_array.push_back(vector<REAL>(n));
    dcsr_mxv(Qj, (vj_array.rbegin()+1)->data(), vj_array.back().data());
  }
  // Finding the optimal basis
  assert(Nj_array.size() == vj_array.size());
  int num_levels = Nj_array.size();
  vector<vector<REAL>> sums(num_levels);
  vector<vector<bool>> labels(num_levels);
  for (int j = num_levels - 1; j >= 0; --j) {
    int i = 0, aux, offset;
    if (j < num_levels - 1) {
      aux = Nj_array[j] - 2 * Nj_array[j+1];
      offset = (1 << (j+1)) * Nj_array[j+1];
    }
    for (int k = 0; k < (1 << j); ++k) {
      sums[j].push_back(0.0);
      for (int l = 0; l < Nj_array[j]; ++l) {
        sums[j].back() += pow(abs(vj_array[j][i++]), p);
      }
      if (j < num_levels - 1) {
        REAL sum_aux = accumulate(
            vj_array[j+1].begin() + offset + k*aux,
            vj_array[j+1].begin() + offset + (k+1)*aux,
            0);
        REAL sum_children = sums[j+1][2*k] + sums[j+1][2*k+1] + sum_aux;
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
  vector<REAL> v_e;
  v_e.reserve(n);
  function<void(int, int)> encode = [&] (int j, int k) {
    if (labels[j][k]) {
      v_e.insert(
          v_e.end(),
          vj_array[j].begin() + k*Nj_array[j],
          vj_array[j].begin() + (k+1)*Nj_array[j]);
    } else {
      encode(j+1, 2*k);
      encode(j+1, 2*k+1);
      int offset = (1 << (j+1)) * Nj_array[j+1];
      int aux = Nj_array[j] - 2 * Nj_array[j+1];
      v_e.insert(
          v_e.end(),
          vj_array[j+1].begin() + offset + k*aux,
          vj_array[j+1].begin() + offset + (k+1)*aux);
    }
  };
  encode(0, 0);
  // Trunk
  for (int j = 0; j < n; ++j) {
    vt_sort[j] = abs(v_e[j]);
  }
  sort(vt_sort.begin(), vt_sort.end(), greater<double>());
  for (int j = 0; j < n; ++j) {
    if (abs(v_e[j]) < vt_sort[threshold]) {
      v_e[j] = 0;
    }
  }
  // Decompress
  auto it = v_e.begin();
  function<vector<REAL>(int, int)> decode = [&] (int j, int k) -> vector<REAL> {
    if (labels[j][k]) {
      auto res = vector<REAL>(it, it + Nj_array[j]);
      it += Nj_array[j];
      return res;
    } else {
      auto v1 = decode(j+1, 2*k);
      auto v2 = decode(j+1, 2*k+1);
      vector<REAL> v_segment;
      int Nj = Nj_array[j];
      v_segment.reserve(Nj);
      v_segment.insert(v_segment.end(), v1.begin(), v1.end());
      v_segment.insert(v_segment.end(), v2.begin(), v2.end());
      int aux = Nj - 2 * Nj_array[j+1];
      v_segment.insert(v_segment.end(), it, it + aux);
      it += aux;
      vector<int> Is, Js;
      for (int i = 0; i < 2*Nj_array[j+1]; ++i) {
        Is.push_back(i);
        Js.push_back(i);
      }
      int offset = (1 << (j+1)) * Nj_array[j+1];
      for (int i = 0; i < aux; ++i) {
        Is.push_back(offset + i);
        Js.push_back(2*Nj_array[j+1] + i);
      }
      dCSRmat *Qj_block = (dCSRmat *)malloc(sizeof(dCSRmat));
      dcsr_getblk(Qj_array[j], Is.data(), Js.data(), Nj, Nj, Qj_block);
      dCSRmat *Qj_block_t = (dCSRmat *)malloc(sizeof(dCSRmat));
      dcsr_trans(Qj_block, Qj_block_t);
      vector<REAL> v_res(Nj);
      dcsr_mxv(Qj_block_t, v_segment.data(), v_res.data());
      return v_res;
    }
  };
  auto v3_vector = decode(0, 0);
  copy(v3_vector.begin(), v3_vector.end(), v3);
  // Compute the error
  REAL e3[n];
  array_axpyz(n, -1.0, v, v3, e3);
  cout << endl << "Adaptive Encoding" << endl
       << "Norm of vector ||v||: " << array_norm2(n, v) << endl
       << "Norm of error  ||v-v3||: " << array_norm2(n, e3) << endl
       << "Relative error ||v-v3||/||v||: "
       << array_norm2(n, e3) / array_norm2(n, v) << endl;
}

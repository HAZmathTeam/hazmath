#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iterator>
#include <numeric>
#include <algorithm>
#include <cassert>
#include <functional>

#include "graph.hpp"

using namespace std;

extern "C" {
  void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda, double *w,
              double *work, int *lwork, int *info);
}

void setup_hierarchy(const char *file, dCSRmat *&A, vector<dCSRmat *> &Qj_array,
    vector<int> &Nj_array);

void comp_decomp(int argc, char *argv[], double *v, dCSRmat *A,
    const vector<dCSRmat *> &Qj_array, const vector<int> &Nj_array,
    double* v2, double *v3);

int main(int argc, char *argv[]) {
  assert(argc > 1);

  int side = 512, nnz = (side-1)*side*2, edge_count = 0, n = side*side;
  ofstream tempfile("temp");
  tempfile << side*side << ' ' << side*side << ' ' << nnz << '\n';
  for (int i = 0; i < side; ++i) {
    for (int j = 0; j < side; ++j) {
      int index = i*side + j + 1;
      if (j < side - 1) {
        tempfile << index << ' ' << index + 1 << '\n';
        ++edge_count;
      }
      if (i < side - 1) {
        tempfile << index << ' ' << index + side << '\n';
        ++edge_count;
      }
    }
  }
  assert(edge_count == nnz);
  tempfile.close();

  dCSRmat *A;
  vector<dCSRmat *> Qj_array;
  vector<int> Nj_array;
  setup_hierarchy("temp", A, Qj_array, Nj_array);

  // Compress/decompress a smooth vector and compute the error
  const string prefix(argv[1]);

  vector<string> months{"01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"};
  for (auto mon : months) {
    string filename = prefix+mon;
    ifstream file(filename.c_str());
    string line;
    vector<vector<double>> matrix;
    while (getline(file, line)) {
      istringstream iss(line);
      matrix.push_back(vector<double>(
          istream_iterator<double>(iss), istream_iterator<double>()));
    }
    file.close();
    for (auto row : matrix) {
      assert(matrix.size() == row.size());
    }

    REAL *v = (REAL *)malloc(sizeof(REAL)*n);
    int v_ind = 0;
    for (auto row : matrix) {
      for (auto entry : row) {
        v[v_ind++] = entry;
      }
    }

    cout << endl << endl << "Month" + mon << endl;
    REAL *v2 = (REAL *)malloc(sizeof(REAL)*n);
    REAL *v3 = (REAL *)malloc(sizeof(REAL)*n);
    comp_decomp(argc-1, argv+1, v, A, Qj_array, Nj_array, v2, v3);

    auto write = [&] (string filename, REAL data[]) {
      ofstream ofs(filename);
      for (int i = 0; i < side; ++i) {
        for (int j = 0; ; ++j) {
          ofs << data[i*side+j];
          if (j == side - 1) break;
          ofs << ' ';
        }
        ofs << '\n';
      }
      ofs.close();
    };
    write(filename+".1", v2);
    write(filename+".2", v3);

    free(v);
    free(v2);
    free(v3);
  }

  dcsr_free(A);
  for (auto Qj : Qj_array) {
    dcsr_free(Qj);
  }
  return 0;
}

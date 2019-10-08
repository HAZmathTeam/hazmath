/* Example to compress precipitation data in Punjab and decompress it,
 * then report the error.
 * Usage:
 *   ./ex2 graphs/Precip/AgMERRA_Precip_by_mo_tot/agp1980_ -k 100
 *   ./ex2 graphs/Precip/AgMERRA_Precip_by_mo_tot/agp1980_ -l
 */
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <cassert>
#include <unistd.h>
#include "graph.h"
#include "algorithm.h"

using namespace std;

int main(int argc, char *argv[]) {
  int threshold = 100;
  bool opt_k = false;
  bool opt_l = false;
  int c;
  opterr = 0;

  while ((c = getopt(argc, argv, "k:l")) != -1) {
    switch (c) {
      case 'k':
        threshold = stoi(optarg);
        opt_k = true;
        break;
      case 'l':
        opt_l = true;
        break;
      case '?':
        if (optopt == 'k') {
          fprintf(stderr, "Option -%c requires an argument.\n", optopt);
        }
        else if (isprint(optopt)) {
          fprintf(stderr, "Unknown option `-%c'.\n", optopt);
        }
        else {
          fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
        }
        return 1;
      default:
        abort();
    }
  }
  if (optind != argc - 1) {
    cerr << "Too many arguments." << endl;
    return 1;
  }
  if (opt_k && opt_l) {
    cerr << "Conflicting options -k and -l." << endl;
    return 1;
  }

  int side = 512, nnz = (side-1)*side*2, edge_count = 0, n = side*side;
  if (threshold > n - 1) {
    threshold = n - 1;
  }
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
  setupHierarchy("temp", A, Qj_array, Nj_array);
  remove("temp");

  // Compress/decompress a smooth vector and compute the error
  const string prefix(argv[optind]);

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
    if (!opt_l) {
      compAndDecomp(v, A, Qj_array, Nj_array, threshold, 1.0, v2, v3);

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
      write(filename+"_1", v2);
      write(filename+"_2", v3);
    }
    else {
      string ofilename = filename;
      ofilename.insert(ofilename.rfind('/'), "/levels");
      ofstream ofs(ofilename + ".data");
      ofs << "# Compression results for plain and adaptive encoding" << endl;
      for (int th = 1; th < n; th <<= 1) {
        compAndDecomp(v, A, Qj_array, Nj_array, th, 1.0, v2, v3);
        double e2[n], e3[n];
        array_axpyz(n, -1.0, v, v2, e2);
        array_axpyz(n, -1.0, v, v3, e3);
        ofs << th << '\t' << array_norm2(n, v) << '\t' << array_norm2(n, e2)
            << '\t' << array_norm2(n, e3) << endl;
      }
      ofs.close();
    }

    free(v2);
    free(v3);
    free(v);
  }

  dcsr_free(A);
  for (auto Qj : Qj_array) {
    dcsr_free(Qj);
  }
  return 0;
}

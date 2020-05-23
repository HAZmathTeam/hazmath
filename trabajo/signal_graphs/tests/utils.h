#include <iterator>
#include <sstream>
#include <string>

template <typename T> std::string arr2str(T *array, int size) {
  std::ostringstream oss;
  oss << '[';
  if (size > 0) {
    std::copy(array, array + size - 1, std::ostream_iterator<T>(oss, ", "));
    oss << array[size - 1];
  }
  oss << ']';
  return oss.str();
}

template <typename T>
void assertArraysEqual(T *first, const std::vector<T> &second,
                       T tol = T(1e-10)) {
  for (int i = 0; i < second.size(); ++i) {
    if (abs(first[i] - second[i]) > tol) {
      throw std::runtime_error('\n' + arr2str(first, second.size()) +
                               " != " + arr2str(second.data(), second.size()));
    }
  }
}

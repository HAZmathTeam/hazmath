#include <iostream>

template<typename T>
void printArray(T* array, int size) {
  std::cout << '[';
  for (int i = 0; i < size; ++i) {
    std::cout << array[i];
    if (i < size - 1) {
      std::cout << ", ";
    }
  }
  std::cout << ']';
}

template<typename T>
void assertArraysEqual(T* first, const std::vector<T>& second) {
  for (int i = 0; i < second.size(); ++i) {
    if (!(first[i] == second[i])) {
      printArray(first, second.size());
      std::cout << " != ";
      printArray(second.data(), second.size());
      std::cout << std::endl;
      throw;
    }
  }
}

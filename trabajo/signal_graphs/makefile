CXX = c++
CXXFLAGS = -Wall -Wfatal-errors -std=c++11
INCLUDE = -I../../include
LIBS = -lhazmath -llapack
HAZLIB = ../../lib

HEADERS = algorithm.h graph.h hamiltonian_algorithm.h hazmath_include.h \
		tests/utils.h utils.h

SOURCE_FILES = algorithm.cpp graph.cpp hamiltonian_algorithm.cpp utils.cpp
OBJECT_FILES = $(SOURCE_FILES:.cpp=.o)

TESTS = get_hamiltonian_path_test get_random_prolongation_test \
		get_walsh_basis_test degree_based_matching_test get_eigenvectors_test
TEST_FILES = $(addprefix tests/, $(TESTS))

EXAMPLES = ex1 ex2


.PHONY: all clean

all: $(EXAMPLES)

%.o: %.cpp $(HEADERS)
	+$(CXX) -c $(CXXFLAGS) $(INCLUDE) -o $@ $<

$(EXAMPLES): % : $(OBJECT_FILES) %.cpp
	+$(CXX) $(CXXFLAGS) $(INCLUDE) $@.cpp $(OBJECT_FILES) -o $@ -L$(HAZLIB) $(LIBS)

# Tests
tests/%: $(OBJECT_FILES) tests/%.cpp
	+$(CXX) $(CXXFLAGS) $(INCLUDE) $@.cpp $(OBJECT_FILES) -o $@ -L$(HAZLIB) $(LIBS)
	+$@

test: $(TEST_FILES)

clean:
	+rm -rf $(EXAMPLES) $(OBJECT_FILES) $(TEST_FILES)

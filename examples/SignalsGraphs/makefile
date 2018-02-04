CXX = mpicxx

SOURCE_FILES = graph.cpp utils.cpp algorithm.cpp
OBJECT_FILES = $(SOURCE_FILES:.cpp=.o)
EXAMPLES = ex1 ex2

HEADERS = graph.hpp
INCLUDE = -I../../include
LIBS = -lhazmath -llapack
HAZLIB = ../../lib

CXXFLAGS = -Wall -Wfatal-errors -std=c++11


.PHONY: all clean

all: $(EXAMPLES)

%.o: %.cpp $(HEADERS)
	+$(CXX) -c $(CXXFLAGS) $(INCLUDE) -o $@ $<

$(EXAMPLES) : % : $(OBJECT_FILES) %.cpp
	+$(CXX) $(CXXFLAGS) $(INCLUDE) $@.cpp $(OBJECT_FILES) -o $@ -L$(HAZLIB) $(LIBS)

clean:
	+rm -rf $(EXAMPLES) $(OBJECT_FILES)

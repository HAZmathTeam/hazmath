CXX = mpicxx

SOURCE_FILES = graph.cpp utils.cpp ex2.cpp
OBJECT_FILES = $(SOURCE_FILES:.cpp=.o)
EXAMPLES = ex2

HEADERS = graph.hpp
INCLUDE = -I../../include -I../../../boost_1_65_1
LIBS = -lhazmath -llapack
HAZLIB = ../../lib

CXXFLAGS = -Wall $(INCLUDE) -std=c++11


.PHONY: all examples clean

all: $(EXAMPLES)

$(EXAMPLES): $(OBJECT_FILES)
	+$(CXX) -L$(HAZLIB) $(OBJECT_FILES) -o ex2 $(LIBS) -Wfatal-errors

%.o: %.cpp $(HEADERS)
	+$(CXX) -c $(CXXFLAGS) -o $@ $< -Wfatal-errors

clean:
	+rm -rf *.o $(EXAMPLES)

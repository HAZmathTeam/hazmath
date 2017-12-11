CXX = g++

SOURCE_FILES = graph.cpp utils.cpp ex4.cpp
OBJECT_FILES = $(SOURCE_FILES:.cpp=.o)
EXAMPLES = ex4

HEADERS = graph.hpp
INCLUDE = ../../include
LIBS = -lhazmath -llapack
HAZLIB = ../../lib

CXXFLAGS = -g -Wall -I$(INCLUDE)


.PHONY: all examples clean

all: $(EXAMPLES)

$(EXAMPLES): $(OBJECT_FILES)
	+$(CXX) -L$(HAZLIB) $(OBJECT_FILES) -o ex4 $(LIBS) -Wfatal-errors

%.o: %.cpp $(HEADERS)
	+$(CXX) -c $(CXXFLAGS) -o $@ $< -Wfatal-errors

clean:
	+rm -rf *.o $(EXAMPLES)

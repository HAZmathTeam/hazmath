CXX = g++
CXXFLAGS = -g -Wall

SOURCE_FILES = graph.cpp
OBJECT_FILES = $(SOURCE_FILES:.cpp=.o)
EXAMPLES = ex4

INCFLAGS = -I../../include
LIBS = -L../../lib -lhazmath -llapacke

.PHONY: all lib examples clean

%.o: %.cpp %.hpp
	$(CXX) -c $(CXXFLAGS) $(INCFLAGS) -o $@ $< -Wfatal-errors

%: %.cpp libgraph.a
	$(CXX) $(CXXFLAGS) $(INCFLAGS) $@.cpp -o $@ -L. -lgraph $(LIBS) -Wfatal-errors

all: lib examples

lib: libgraph.a

libgraph.a: $(OBJECT_FILES)
	ar csr libgraph.a $(OBJECT_FILES)
	ranlib libgraph.a

examples: $(EXAMPLES)

clean:
	rm libgraph.a $(OBJECT_FILES) $(EXAMPLES)

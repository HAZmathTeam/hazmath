CXX = g++
CXXFLAGS = -g -Wall

SOURCE_FILES = graph.cpp
OBJECT_FILES = $(SOURCE_FILES:.cpp=.o)
OBJ22 = ex4.o eig_call_fortran.o
EXAMPLES = ex4

INCFLAGS = -I../../include
LIBS = -L../../lib  -llapack -lhazmath

.PHONY: all lib examples clean

%.o: %.cpp %.hpp
	$(CXX) -c $(CXXFLAGS) $(INCFLAGS) -o $@ $<  -Wfatal-errors

$(EXAMPLES): $(OBJ22)	libgraph.a
	$(CXX) $(CXXFLAGS) $(INCFLAGS)  $(OBJ22) -o ex4 -L. -lgraph $(LIBS) -Wfatal-errors

all: lib examples

lib: libgraph.a

libgraph.a: $(OBJECT_FILES)
	ar csr libgraph.a $(OBJECT_FILES)
	ranlib libgraph.a

examples: $(EXAMPLES)

clean:
	+rm -rf libgraph.a $(OBJECT_FILES) $(EXAMPLES)

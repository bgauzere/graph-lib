IDIR = ./include
LSAPE_DIR=/home/bgauzere/dev/lsape/include/
EIGEN_DIR=/usr/include/eigen3/
CXXFLAGS = -I$(IDIR) -I$(LSAPE_DIR) -I$(EIGEN_DIR)  -Wall  -std=c++11
BINDIR = ./bin
TESTDIR = ./test


ODIR = ./obj
SRCDIR = ./src
_DEPS = graph.h  utils.h SymbolicGraph.h GraphEditDistance.h ConstantGraphEditDistance.h Dataset.h BipartiteGraphEditDistance.h RandomWalksGraphEditDistance.h IPFPGraphEditDistance.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_SRCDEPS = utils.cpp  SymbolicGraph.cpp ConstantGraphEditDistance.cpp RandomWalksGraphEditDistance.cpp
DEPS_SRC += $(patsubst %,$(SRCDIR)/%,$(_DEPS_SRC))

_OBJ = utils.o SymbolicGraph.o ConstantGraphEditDistance.o RandomWalksGraphEditDistance.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

# all: $(BINDIR)/test_GraphEditDistance $(BINDIR)/contestGraphEditDistance
all:$(TESTDIR)/test_graph

debug: CXXFLAGS += -DDEBUG -g
debug: all

optim: CXXFLAGS += -O3
optim: all

$(TESTDIR)/test_graph: $(DEPS) $(OBJ) $(TESTDIR)/test_graph.cpp
	$(CXX) -o $@ $^ $(CXXFLAGS) -ltinyxml

$(BINDIR)/%: $(OBJ) $(SRCDIR)/%.cpp
	$(CXX) -o $@ $^ $(CXXFLAGS)

$(ODIR)/%.o: $(SRCDIR)/%.cpp $(DEPS) $(DEPS_SRC)
	$(CXX) -c -o $@ $< $(CXXFLAGS)


.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~  $(BINDIR)/*

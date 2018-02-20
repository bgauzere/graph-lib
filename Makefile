IDIR = ./include

## LSAPE directory 
#LSAPE_DIR=$(LSAPE_DIR)#/home/bgauzere/Téléchargements/lsape/include/

## Eigen library
EIGEN_DIR=/usr/include/eigen3/

CXXFLAGS = -I$(IDIR) -I$(LSAPE_DIR) -I$(EIGEN_DIR) -Wall  -std=c++11 -g #-Werror

BINDIR = ./bin
TESTDIR = ./test
ODIR = ./obj
SRCDIR = ./src

_DEPS = graph.h  utils.h SymbolicGraph.h GraphEditDistance.h ConstantGraphEditDistance.h Dataset.h MultiGed.h BipartiteGraphEditDistance.h BipartiteGraphEditDistanceMulti.h RandomWalksGraphEditDistance.h RandomWalksGraphEditDistanceMulti.h IPFPGraphEditDistance.h  MultistartRefinementGraphEditDistance.h IPFPZetaGraphEditDistance.h  GNCCPGraphEditDistance.h CMUCostFunction.h CMUGraph.h  CMUDataset.h LetterCostFunction.h LetterGraph.h LetterDataset.h BipartiteLowerBoundGraphEditDistance.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_SRCDEPS = utils.cpp  SymbolicGraph.cpp ConstantGraphEditDistance.cpp RandomWalksGraphEditDistance.cpp RandomWalksGraphEditDistanceMulti.cpp  CMUCostFunction.cpp  CMUGraph.cpp  CMUDataset.cpp  LetterCostFunction.cpp  LetterGraph.cpp LetterDataset.cpp 
DEPS_SRC += $(patsubst %,$(SRCDIR)/%,$(_DEPS_SRC))

_OBJ = utils.o SymbolicGraph.o ConstantGraphEditDistance.o RandomWalksGraphEditDistance.o CMUCostFunction.o CMUGraph.o CMUDataset.o LetterGraph.o LetterCostFunction.o LetterDataset.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

_OBJ_QAP = utils.o QAPLibGraph.o QAPLibCostFunction.o QAPLibDataset.o
OBJ_QAP = $(patsubst %,$(ODIR)/%,$(_OBJ_QAP))

# all: $(BINDIR)/test_GraphEditDistance $(BINDIR)/contestGraphEditDistance
all:$(TESTDIR)/test_graph $(TESTDIR)/chemical-edit-distances #$(TESTDIR)/benchmark

debug: CXXFLAGS += -DDEBUG -g
debug: $(TESTDIR)/chemical-lower-bounds

with_times: CXXFLAGS += -D PRINT_TIMES
with_times: $(TESTDIR)/chemical-edit-distances

multithread: CXXFLAGS += -fopenmp
multithread: $(TESTDIR)/chemical-edit-distances

multithread_with_times: CXXFLAGS += -fopenmp -D PRINT_TIMES
multithread_with_times: $(TESTDIR)/chemical-edit-distances

qaplib_mt: CXXFLAGS += -fopenmp
qaplib_mt: $(TESTDIR)/QAPlib

optim: CXXFLAGS += -O3
#optim: all

$(TESTDIR)/xp_PRL2017: $(DEPS) $(OBJ) $(TESTDIR)/xp_PRL2017.cpp
	$(CXX) -o $@ $^ $(CXXFLAGS) -ltinyxml

# $(TESTDIR)/benchmark: $(DEPS) $(OBJ) $(TESTDIR)/benchmark.cpp
# 	$(CXX) -o $@ $^ $(CXXFLAGS) -ltinyxml

$(TESTDIR)/chemical-edit-distances: $(TESTDIR)/computeDistances.cpp $(OBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS) -ltinyxml

$(TESTDIR)/chemical-lower-bounds: $(TESTDIR)/computeLowerBounds.cpp $(OBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS) -ltinyxml

$(TESTDIR)/QAPlib: $(TESTDIR)/QAPLib.cpp $(OBJ_QAP)
	$(CXX)  -D PRINT_TIMES -o $@ $^ $(CXXFLAGS) -ltinyxml

$(TESTDIR)/test_graph: $(DEPS) $(OBJ) $(TESTDIR)/test_graph.cpp
	$(CXX) -o $@ $^ $(CXXFLAGS) -ltinyxml

$(BINDIR)/%: $(OBJ) $(SRCDIR)/%.cpp
	$(CXX) -o $@ $^ $(CXXFLAGS)

$(ODIR)/%.o: $(SRCDIR)/%.cpp $(DEPS) $(DEPS_SRC)
	$(CXX) -c -o $@ $< $(CXXFLAGS)


.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~  $(BINDIR)/*

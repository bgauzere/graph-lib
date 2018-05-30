IDIR = ../include

## LSAPE directory
#  LSAPE_DIR=/export/home/daller/Documents/dev/lsape/include/

## Eigen library
EIGEN_DIR=/usr/include/eigen3/

## Graph-lib static library path
STATIC_GRAPH_LIB=../bin/graphlib.a

CXXFLAGS = -I$(IDIR) -I$(LSAPE_DIR) -I$(EIGEN_DIR) -Wall  -std=c++11 -O3 -fopenmp

TESTDIR = ./
ODIR = ../obj
SRCDIR = ../src

_OBJ = gl_utils.o SymbolicGraph.o ConstantGraphEditDistance.o RandomWalksGraphEditDistance.o CMUCostFunction.o CMUGraph.o CMUDataset.o LetterGraph.o LetterCostFunction.o LetterDataset.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


all: $(TESTDIR)/frontend


$(TESTDIR)/frontend_example: $(TESTDIR)/frontend_example.cpp $(STATIC_GRAPH_LIB)
	$(CXX) -o $@ $^ $(CXXFLAGS) -ltinyxml

$(STATIC_GRAPH_LIB): $(OBJ)
	ar rcs $@ $^

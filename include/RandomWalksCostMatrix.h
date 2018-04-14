/**
 * @file RandomWalksCostMatrix.h
 * @author Benoit <<benoit.gauzere@insa-rouen.fr>> 
 * @version     0.0.1 - Wed Feb  8 2017
 *
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *
 * Description of the program objectives.
 * All necessary references.
 */

#ifndef __RANDOMWALKSGRAPHEDITDISTANCE_H__
#define __RANDOMWALKSGRAPHEDITDISTANCE_H__

#include <Eigen/Dense>
using namespace Eigen;

#include "RiesenCostMatrix.h"
#include "ConstantGraphEditDistance.h"

#include "SymbolicGraph.h"

class RandomWalksCostMatrix:
public virtual RiesenCostMatrix<int, int>
{
protected:
  ConstantEditDistanceCost * _cf;
  int _k;

  static int * labeledKron(int *m1, int nb_rows_m1,int nb_cols_m1,
			   int * m2, int nb_rows_m2, int nb_cols_m2,
			   int sizeWx[2]);
  static MatrixXi histoLab(int nbLab,  RowVectorXi IL, MatrixXi W);

  void computeCostMatrix(Graph<int,int> * g1,
			 Graph<int,int> * g2);

public:
  RandomWalksCostMatrix(ConstantEditDistanceCost * costFunction, int k):
    RiesenCostMatrix<int,int>(costFunction),_cf(costFunction),_k(k){};
};

#endif // __RANDOMWALKSCOSTMATRIX_H__

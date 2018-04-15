/**
 * @file CostMatrix.h
 * @author bgauzere <<benoit.gauzere@insa-rouen.fr>> 
 * @version     0.0.1 - Fri Apr 13 2018
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 */

#ifndef __COSTMATRIX_H__
#define __COSTMATRIX_H__
#include "graph.h"
#include "GraphEditDistance.h"

template<class NodeAttribute, class EdgeAttribute>
class CostMatrix
{
protected:
  EditDistanceCost<NodeAttribute,EdgeAttribute> * _cf;
  double * _C;
  int _n; //number of rows of C
  int _m; //number of columnsof C
  
public:
  CostMatrix(EditDistanceCost<NodeAttribute,EdgeAttribute> * cf):_cf(cf),_C(NULL),_n(0),_m(0){};
  
  virtual void computeCostMatrix(Graph<NodeAttribute,EdgeAttribute> * g1,
				 Graph<NodeAttribute,EdgeAttribute> * g2){};
  virtual double * getC(){return this->_C;}

  ~CostMatrix(){
    if (this->_C)
      delete [] _C;
  }

};

#endif // __COSTMATRIX_H__

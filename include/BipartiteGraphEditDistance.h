/**
 * @file BipartiteGraphEditDistance.h
 * @author Benoit <<benoit.gauzere@insa-rouen.fr>> 
 * @version     0.0.1 - Mon Feb  6 2017
 *
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *
 * Bipartite Graph edit distance algorithm as described in Riesen's
 * book.  Riesen, K. (2015). Structural pattern recognition with graph
 * edit distance. Advances in Computer Vision and Pattern Recognition,
 * Cham.
 */

#ifndef __BIPARTITEGRAPHEDITDISTANCE_H__
#define __BIPARTITEGRAPHEDITDISTANCE_H__

#include "GraphEditDistance.h"
#include "gl_utils.h"
#include "CostMatrix.h"
#include "LinearSolver.h"

/*TODO : 
   - donner la possibilité de récupérer le mapping ?
   - Spécifier le solver
*/

template<class NodeAttribute, class EdgeAttribute>
class BipartiteGraphEditDistance:
  public GraphEditDistance<NodeAttribute, EdgeAttribute>
{
private:

protected:
  LinearSolver<NodeAttribute, EdgeAttribute> * _solver;
  CostMatrix<NodeAttribute,EdgeAttribute> * _CM;
  

public:
  BipartiteGraphEditDistance(EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction,
			     CostMatrix<NodeAttribute,EdgeAttribute> * CM,
			     LinearSolver<NodeAttribute, EdgeAttribute> * solver):
    GraphEditDistance<NodeAttribute,EdgeAttribute>(costFunction),_solver(solver),_CM(CM)
  {};

  double getTimedOptimalMapping(Graph<NodeAttribute,EdgeAttribute> * g1,
				Graph<NodeAttribute,EdgeAttribute> * g2,
				unsigned int * G1_to_G2, unsigned int * G2_to_G1);
  virtual void getOptimalMapping(Graph<NodeAttribute,EdgeAttribute> * g1,
				 Graph<NodeAttribute,EdgeAttribute> * g2,
				 unsigned int * G1_to_G2, unsigned int * G2_to_G2);
  double getUpperBound(Graph<NodeAttribute,EdgeAttribute> * g1,
		       Graph<NodeAttribute,EdgeAttribute> * g2,
		       double & time);
  

  virtual ~BipartiteGraphEditDistance(){
    ;    
  }

  virtual BipartiteGraphEditDistance<NodeAttribute, EdgeAttribute> * clone() const
  {
    return new BipartiteGraphEditDistance<NodeAttribute, EdgeAttribute> (*this);
  }

};


//A remonter a GraphEditDistance.h ?
template<class NodeAttribute, class EdgeAttribute>
double BipartiteGraphEditDistance<NodeAttribute, EdgeAttribute>::
getTimedOptimalMapping(Graph<NodeAttribute,EdgeAttribute> * g1,
		       Graph<NodeAttribute,EdgeAttribute> * g2,
		       unsigned int * G1_to_G2, unsigned int * G2_to_G1){
  int n=g1->Size();
  int m=g2->Size();
  // Compute C
  this->_CM->computeCostMatrix(g1,g2);
  //Compute optimal assignement
  clock_t t = clock();
  this->_solver->solve(this->_CM->getC(),n+1,m+1, g1,g2,G1_to_G2, G2_to_G1);
  t = clock() - t;
  return ((float)t) / CLOCKS_PER_SEC;

}

template<class NodeAttribute, class EdgeAttribute>
void BipartiteGraphEditDistance<NodeAttribute, EdgeAttribute>::
getOptimalMapping(Graph<NodeAttribute,EdgeAttribute> * g1,
		  Graph<NodeAttribute,EdgeAttribute> * g2,
		  unsigned int * G1_to_G2, unsigned int * G2_to_G1){
  int n=g1->Size();
  int m=g2->Size();
  // Compute C
  this->_CM->computeCostMatrix(g1,g2);
  this->_solver->solve(this->_CM->getC(),n+1, m+1, g1,g2,G1_to_G2, G2_to_G1);
}

template<class NodeAttribute, class EdgeAttribute>
double BipartiteGraphEditDistance<NodeAttribute, EdgeAttribute>::
getUpperBound(Graph<NodeAttribute,EdgeAttribute> * g1,
	      Graph<NodeAttribute,EdgeAttribute> * g2,double & comptime){
  int n=g1->Size();
  int m=g2->Size();
  unsigned int * G1_to_G2 = new unsigned int[n+1];
  unsigned int * G2_to_G1 = new unsigned int[m+1];
  comptime = this->getTimedOptimalMapping(g1,g2,G1_to_G2,G2_to_G1);
  double ged = this->GedFromMapping(g1,g2,G1_to_G2,n,G2_to_G1,m);
  delete [] G1_to_G2;
  delete [] G2_to_G1;
  return ged;
}


#endif // __BIPARTITEGRAPHEDITDISTANCE_H__

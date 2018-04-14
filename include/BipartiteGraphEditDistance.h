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

#include "lsape.h"
#include "GraphEditDistance.h"
#include "gl_utils.h"
#include "CostMatrix.h"


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
  enum lsape::LSAPE_MODEL  my_solver;
  CostMatrix<NodeAttribute,EdgeAttribute> * _CM;


public:
  BipartiteGraphEditDistance(EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction,
			     CostMatrix<NodeAttribute,EdgeAttribute> * CM):
    GraphEditDistance<NodeAttribute,EdgeAttribute>(costFunction),my_solver(lsape::ECBP),_CM(CM)
  {};

  BipartiteGraphEditDistance(EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction,
			     CostMatrix<NodeAttribute,EdgeAttribute> * CM,
			     enum lsape::LSAPE_MODEL f):
    GraphEditDistance<NodeAttribute,EdgeAttribute>(costFunction),
    my_solver(f),_CM(CM){};

  void setSolver(enum lsape::LSAPE_MODEL f){
    this->my_solver = f;
  }

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
  double *u = new double[n+1];
  double *v = new double[m+1];
  clock_t t = clock();
  lsape::lsapeSolver<double>(this->_CM->getC(),n+1,m+1, G1_to_G2, G2_to_G1, u,v,my_solver,1);
  t = clock() - t;
  delete [] u;
  delete [] v;
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
  //Compute optimal assignement
  double *u = new double[n+1];
  double *v = new double[m+1];
  lsape::lsapeSolver<double>(this->_CM->getC(),n+1,m+1, G1_to_G2, G2_to_G1, u,v,my_solver,1);
  delete [] u;
  delete [] v;

}

template<class NodeAttribute, class EdgeAttribute>
double BipartiteGraphEditDistance<NodeAttribute, EdgeAttribute>::
getUpperBound(Graph<NodeAttribute,EdgeAttribute> * g1,
	      Graph<NodeAttribute,EdgeAttribute> * g2,double & comptime){
  int n=g1->Size();
  int m=g2->Size();
  unsigned int * G1_to_G2 = new unsigned int[n];
  unsigned int * G2_to_G1 = new unsigned int[m];
  comptime = this->getTimedOptimalMapping(g1,g2,G1_to_G2,G2_to_G1);
  double ged = this->GedFromMapping(g1,g2,G1_to_G2,n,G2_to_G1,m);
  delete [] G1_to_G2;
  delete [] G2_to_G1;
  return ged;
}


#endif // __BIPARTITEGRAPHEDITDISTANCE_H__

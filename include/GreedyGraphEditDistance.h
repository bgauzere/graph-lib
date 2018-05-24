/**
 * @file GreedyBipartiteGraphEditDistance.h
 * @author Evariste <<evariste.daller@unicaen.fr>>
 *
 *
 */


#ifndef __GREEDYBIPARTITEGED_H__
#define __GREEDYBIPARTITEGED_H__


#include "BipartiteGraphEditDistanceMulti.h"


template <class NodeAttribute, class EdgeAttribute>
class GreedyGraphEditDistance :
  public BipartiteGraphEditDistanceMulti<NodeAttribute, EdgeAttribute>
{

protected:

  //int* _costs;

public:

  GreedyGraphEditDistance( EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction,
			   int nep ):
    BipartiteGraphEditDistance<NodeAttribute,EdgeAttribute>(costFunction),
    BipartiteGraphEditDistanceMulti<NodeAttribute,EdgeAttribute>(costFunction, nep)
  {};


public:

  /**
   * @brief Allocate and retruns `k` optimal mappings between `g1` and `g2`
   * @param k  The number of mappings to compute, -1 to get all perfect matchings
   * @param C  The cost matrix
   * @return  A list of mappings given as arrays of int. For each mapping M, `M[i]` is the mapping, in g2, of node i in g1
   * @note  Each array is allocated here and have to be deleted manually
   */
  virtual std::list<unsigned int*> getMappings(Graph<NodeAttribute,EdgeAttribute> * g1,
				      Graph<NodeAttribute,EdgeAttribute> * g2,
				      int k );

};






template<class NodeAttribute, class EdgeAttribute>
std::list<unsigned int*> GreedyGraphEditDistance<NodeAttribute, EdgeAttribute>::
getMappings( Graph<NodeAttribute,EdgeAttribute> * g1,
	     Graph<NodeAttribute,EdgeAttribute> * g2,
	     int k)
{
  unsigned int n=g1->Size();
  unsigned int m=g2->Size();
  
  this->computeCostMatrix(g1, g2);
  
  // Compute integer cost matrix
  int* Ci = new int[(n+1)*(m+1)];

  for (unsigned int j=0; j<=m; j++){
    for (unsigned int i=0; i<=n; i++){
      Ci[sub2ind(i,j,n+1)] = (int)(this->C[sub2ind(i,j,n+1)]);
    }
  }
  // the returned mappings
  std::list<unsigned int*> mappings;

  lsape::cDigraph<unsigned int> dg = lsape::greedySortDigraph<int, unsigned int>(Ci, n+1, m+1);
  lsape::AllPerfectMatchingsEC<unsigned int> apm(dg, n, m, mappings);
  apm.enumPerfectMatchings(dg, this->_nep);
  //mappings = apm.getPerfectMatchings();
  
  delete [] Ci;

  return mappings;
}


#endif

/**
 * @file BipartiteGraphEditDistanceMulti.h
 * @author Évariste <<evariste.daller@unicaen.fr>>
 * @version     0.0.1 - Mon Mar  27 2017
 *
 */

#ifndef __BIPARTITEGRAPHEDITDISTANCEMULTI_H__
#define __BIPARTITEGRAPHEDITDISTANCEMULTI_H__


#include "BipartiteGraphEditDistance.h"
#include "MultiGed.h"
#include "hungarian-lsap.hh"


/**
 * @brief Bipartite version of MultiGed based on the Bunke method for linear cost matrix generation
 */
template<class NodeAttribute, class EdgeAttribute>
class BipartiteGraphEditDistanceMulti :
      public virtual BipartiteGraphEditDistance<NodeAttribute, EdgeAttribute>,
      public MultiGed<NodeAttribute, EdgeAttribute>
{

public:

  BipartiteGraphEditDistanceMulti( EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction,
                                   int _k ):
    BipartiteGraphEditDistance<NodeAttribute,EdgeAttribute>(costFunction),
    MultiGed<NodeAttribute,EdgeAttribute>(_k)
  {};


  BipartiteGraphEditDistanceMulti (
       const BipartiteGraphEditDistanceMulti<NodeAttribute,EdgeAttribute> & other
       ):
    BipartiteGraphEditDistance<NodeAttribute,EdgeAttribute>(other.cf),
    MultiGed<NodeAttribute,EdgeAttribute>(other._nep)
  {
    this->_ged = other._ged;
    this->C = NULL;
  }


public:

  /**
   * @brief compute an optimal mapping between `g1` and `g2`
   *        from k different optimal mappings by minimizing the ged optained
   * @note The GED is computed and set in `ged`
   */
  virtual void getOptimalMapping( Graph<NodeAttribute,EdgeAttribute> * g1,
                                  Graph<NodeAttribute,EdgeAttribute> * g2,
                                  int * G1_to_G2, int * G2_to_G1 );


  virtual std::list<int*> getMappings( Graph<NodeAttribute,EdgeAttribute> * g1,
                                       Graph<NodeAttribute,EdgeAttribute> * g2,
                                       int k = -1 );
  /**
   * @brief Compute the Graph Edit Distance between `g1` and `g2` considering $k$ edit paths
   * @param k  The number of edit paths to compute
   */
  virtual double operator() (Graph<NodeAttribute,EdgeAttribute> * g1,
                              Graph<NodeAttribute,EdgeAttribute> * g2,
                             const int& k=-1);

  /**
   * @brief Compute the GED between `g1` and `g2` as the minimum GED found trough all edit paths
   * @return  calls to operator() (g1, g2, k=1)
   */
  virtual double operator() (Graph<NodeAttribute,EdgeAttribute> * g1,
                             Graph<NodeAttribute,EdgeAttribute> * g2);


  virtual BipartiteGraphEditDistanceMulti<NodeAttribute, EdgeAttribute>* clone() const {
    return new BipartiteGraphEditDistanceMulti<NodeAttribute, EdgeAttribute>(*this);
  }

  virtual ~BipartiteGraphEditDistanceMulti(){
    if (this->C != NULL) {
      delete [] this->C;
      this->C = NULL;
    }
  }
};


//---



template<class NodeAttribute, class EdgeAttribute>
void BipartiteGraphEditDistanceMulti<NodeAttribute, EdgeAttribute>::
getOptimalMapping (Graph<NodeAttribute,EdgeAttribute> * g1,
                   Graph<NodeAttribute,EdgeAttribute> * g2,
                   int * G1_to_G2, int * G2_to_G1 )
{
  delete [] this->C;
  this->computeCostMatrix(g1, g2);
  this->computeOptimalMapping(this, g1, g2, this->C, G1_to_G2, G2_to_G1);
}


template<class NodeAttribute, class EdgeAttribute>
std::list<int*> BipartiteGraphEditDistanceMulti<NodeAttribute, EdgeAttribute>::
getMappings( Graph<NodeAttribute,EdgeAttribute> * g1,
             Graph<NodeAttribute,EdgeAttribute> * g2,
             int k )
{
  if (k == -1) k = this->_nep;
  this->computeCostMatrix(g1, g2);
  std::list<int*> maps = MultiGed<NodeAttribute, EdgeAttribute>::getKOptimalMappings(this->C, g1->Size(), g2->Size(), k);

  delete [] this->C; this->C = NULL;
  return maps;
}


template<class NodeAttribute, class EdgeAttribute>
double BipartiteGraphEditDistanceMulti<NodeAttribute, EdgeAttribute>::
operator() (Graph<NodeAttribute,EdgeAttribute> * g1,
            Graph<NodeAttribute,EdgeAttribute> * g2,
            const int& k)
{
  int n=g1->Size();
  int m=g2->Size();

  int* G1_to_G2 = new int[n];
  int* G2_to_G1 = new int[m];

  if (this->_nep != k) this->_nep = k;
  getOptimalMapping(g1, g2, G1_to_G2, G2_to_G1);

  delete [] G2_to_G1;
  delete [] G1_to_G2;

  return this->_ged;
}


template<class NodeAttribute, class EdgeAttribute>
double BipartiteGraphEditDistanceMulti<NodeAttribute, EdgeAttribute>::
operator() (Graph<NodeAttribute,EdgeAttribute> * g1,
            Graph<NodeAttribute,EdgeAttribute> * g2)
{
  return operator() (g1,g2,this->_nep);
}

#endif

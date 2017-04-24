/**
 * @file BipartiteGraphEditDistanceMulti.h
 * @author Évariste <<evariste.daller@unicaen.fr>>
 * @version     0.0.1 - Mon Mar  27 2017
 *
 * @bug  some local memory alloc changes other functions behavior, see line 303
 */

#ifndef __BIPARTITEGRAPHEDITDISTANCEMULTI_H__
#define __BIPARTITEGRAPHEDITDISTANCEMULTI_H__


#if XP_OUTPUT
#include <ctime>
#endif

#include "BipartiteGraphEditDistance.h"
#include "MultiGed.h"
#include "AllPerfectMatchings-ec.h"
#include "hungarian-lsap.hh"


template<class NodeAttribute, class EdgeAttribute>
class BipartiteGraphEditDistanceMulti :
      public BipartiteGraphEditDistance<NodeAttribute, EdgeAttribute>,
      public MultiGed<NodeAttribute, EdgeAttribute>
{

public:

  BipartiteGraphEditDistanceMulti( EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction,
                                   int _k ):
    BipartiteGraphEditDistance<NodeAttribute,EdgeAttribute>(costFunction),
    MultiGed<NodeAttribute,EdgeAttribute>(_k)
  {};


public:

  /**
   * @brief compute an optimal mapping between <code>g1</code> and <code>g2</code>
   *        from k different optimal mappings by minimizing the ged optained
   * @note The GED is computed and set in <code>ged</code>
   */
  virtual void getOptimalMapping( Graph<NodeAttribute,EdgeAttribute> * g1,
                                  Graph<NodeAttribute,EdgeAttribute> * g2,
                                  int * G1_to_G2, int * G2_to_G1 );

  /**
   * @brief Compute the Graph Edit Distance between <code>g1</code> and <code>g2</code> considering $k$ edit paths
   * @param k  The number of edit paths to compute
   */
  virtual double operator() (Graph<NodeAttribute,EdgeAttribute> * g1,
                              Graph<NodeAttribute,EdgeAttribute> * g2,
                             const int& k=-1);

  /**
   * @brief Compute the GED between <code>g1</code> and <code>g2</code> as the minimum GED found trough all edit paths
   * @return  calls to operator() (g1, g2, k=1)
   */
  virtual double operator() (Graph<NodeAttribute,EdgeAttribute> * g1,
                             Graph<NodeAttribute,EdgeAttribute> * g2);


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
  int n=g1->Size();
  int m=g2->Size();

  delete [] this->C;     //this->C = NULL;

#if XP_OUTPUT
  clock_t t = clock();
#endif
  this->computeCostMatrix(g1, g2);

#if XP_OUTPUT
  t = clock() - t;
  std::cout << ((float)t) / CLOCKS_PER_SEC << ":";
#endif

  this->computeOptimalMapping(this, g1, g2, this->C, G1_to_G2, G2_to_G1);
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

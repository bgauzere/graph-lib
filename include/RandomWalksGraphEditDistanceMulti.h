/**
 * @file RandomWalksGraphEditDistanceMulti.h
 * @author Évariste <<evariste.daller@unicaen.fr>>
 * @version     0.0.1 - Mon Mar  29 2017
 */

#ifndef __RANDOMWALKSGRAPHEDITDISTANCEMULTI_H__
#define __RANDOMWALKSGRAPHEDITDISTANCEMULTI_H__

#include <list>
#include "RandomWalksGraphEditDistance.h"


class RandomWalksGraphEditDistanceMulti :
      public RandomWalksGraphEditDistance
{
private:
  int nep;
  double* Clsap;
  double ged;

public:

  RandomWalksGraphEditDistanceMulti(ConstantEditDistanceCost * costFunction, int k, int nep):
    RandomWalksGraphEditDistance(costFunction, k),
    nep(nep),
    ged(-1.0)
  {}

protected:

    /**
     * @brief Compute the cost matrix corresponding to the lsap given C
     */
    virtual void computeCostMatrixLSAP(int n, int m);

public:

  double getGED() {return ged; } //!< if the returned value is -1, the ged has not been computed @see getOptimalMapping

  /**
   * @brief compute an optimal mapping between <code>g1</code> and <code>g2</code>
   *        from k different optimal mappings by minimizing the ged optained
   * @note The GED is computed and set in <code>ged</code>
   */
  virtual void getOptimalMapping( Graph<int,int> * g1,
                                  Graph<int,int> * g2,
                                  int * G1_to_G2, int * G2_to_G1 );


  /**
   * @brief Allocate and retruns $k$ optimal mappings between <code>g1</code> and <code>g2</code>
   * @param k  The number of mappings to compute, -1 to get all perfect matchings
   * @return  A list of mappings given as arrays of int. For each mapping M, <code>M[i]</code> is the mapping, in g2, of node i in g1
   * @note  Each array is allocated here and have to be deleted manually
   */
  virtual std::list<int*> getKOptimalMappings(Graph<int,int> * g1,
                                              Graph<int,int> * g2,
                                              const int& k=-1);

  /**
   * @brief Compute the Graph Edit Distance between <code>g1</code> and <code>g2</code> considering $k$ edit paths
   * @param k  The number of edit paths to compute
   */
  virtual double operator() (Graph<int,int> * g1,
   			                     Graph<int,int> * g2,
                             const int& k=-1);

  /**
   * @brief Compute the GED between <code>g1</code> and <code>g2</code> as the minimum GED found trough all edit paths
   * @return  calls to operator() (g1, g2, k=1)
   */
  virtual double operator() (Graph<int,int> * g1,
                             Graph<int,int> * g2);
  //*/

  virtual ~RandomWalksGraphEditDistanceMulti();
};

#endif // __RANDOMWALKSGRAPHEDITDISTANCEMULTI_H__

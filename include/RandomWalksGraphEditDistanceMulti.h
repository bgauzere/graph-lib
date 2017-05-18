/**
 * @file RandomWalksGraphEditDistanceMulti.h
 * @author Évariste <<evariste.daller@unicaen.fr>>
 * @version     0.0.1 - Mon Mar  29 2017
 */

#ifndef __RANDOMWALKSGRAPHEDITDISTANCEMULTI_H__
#define __RANDOMWALKSGRAPHEDITDISTANCEMULTI_H__

#include <list>
#include "RandomWalksGraphEditDistance.h"
#include "MultiGed.h"

// activate output for experiments
//#define XP_OUTPUT 1
// Perform xx time tests to get an average
//#define XP_TIME_SAMPLES 1

#if XP_OUTPUT
#include <ctime>
#endif

class RandomWalksGraphEditDistanceMulti :
      public RandomWalksGraphEditDistance,
      public MultiGed<int,int>
{

public:

  RandomWalksGraphEditDistanceMulti(ConstantEditDistanceCost * costFunction, int k, int nep):
    RandomWalksGraphEditDistance(costFunction, k),
    MultiGed<int,int>(nep)
  {}


public:

  /**
   * @brief compute an optimal mapping between <code>g1</code> and <code>g2</code>
   *        from k different optimal mappings by minimizing the ged optained
   * @note The GED is computed and set in <code>ged</code>
   */
  virtual void getOptimalMapping( Graph<int,int> * g1,
                                  Graph<int,int> * g2,
                                  int * G1_to_G2, int * G2_to_G1 );

  /**
   * @brief Compute the Graph Edit Distance between <code>g1</code> and <code>g2</code> considering $k$ edit paths
   * @param k  The number of edit paths to compute
   */
  virtual double operator() (Graph<int,int> * g1,
                             Graph<int,int> * g2);


  virtual ~RandomWalksGraphEditDistanceMulti();
};

#endif // __RANDOMWALKSGRAPHEDITDISTANCEMULTI_H__

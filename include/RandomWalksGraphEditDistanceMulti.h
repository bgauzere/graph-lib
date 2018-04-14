/**
 * @file RandomWalksGraphEditDistanceMulti.h
 * @author Évariste <<evariste.daller@unicaen.fr>>
 * @version     0.0.1 - Mon Mar  29 2017
 */

#ifndef __RANDOMWALKSGRAPHEDITDISTANCEMULTI_H__
#define __RANDOMWALKSGRAPHEDITDISTANCEMULTI_H__

#include <list>
#include "RandomWalksGraphEditDistance.h"
#include "BipartiteGraphEditDistanceMulti.h"
#include "BipartiteGraphEditDistance.h"

// activate output for experiments
//#define XP_OUTPUT 1
// Perform xx time tests to get an average
//#define XP_TIME_SAMPLES 1

#if XP_OUTPUT
#include <ctime>
#endif

/**
 * @brief Multiple solution version of RandomWalksGraphEditDistance
 */
class RandomWalksGraphEditDistanceMulti :
  public BipartiteGraphEditDistanceMulti<int,int>,
  public RandomWalksGraphEditDistance{
  
public:

  RandomWalksGraphEditDistanceMulti(ConstantEditDistanceCost * costFunction, int k, int nep):
    BipartiteGraphEditDistance<int,int>(costFunction),
    BipartiteGraphEditDistanceMulti<int,int>(costFunction, nep),
    RandomWalksGraphEditDistance(costFunction, k)

  {};
};

#endif // __RANDOMWALKSGRAPHEDITDISTANCEMULTI_H__

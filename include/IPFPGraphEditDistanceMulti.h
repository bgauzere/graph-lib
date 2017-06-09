/**
 * @file IPFPGraphEditDistanceMulti.h
 * @author Evariste <<evariste.daller@unicaen.fr>>
 * @version     Jun 9 2017
 *
 */

#ifndef __IPFPGRAPHEDITDISTANCEMULTI_H__
#define __IPFPGRAPHEDITDISTANCEMULTI_H__

#include <list>
#include "IPFPGraphEditDistance.h"
#include "BipartiteGraphEditDistanceMulti.h"
#include "MultiGed.h"


template<class NodeAttribute, class EdgeAttribute>
class IPFPGraphEditDistanceMulti:
  public virtual IPFPGraphEditDistance<NodeAttribute, EdgeAttribute>,
  public MultiGed<NodeAttribute, EdgeAttribute>
{

protected:

  BipartiteGraphEditDistanceMulti<NodeAttribute,EdgeAttribute> * _ed_multi;  //! The multiple solution based GED used (ie. Bunke or Random Walks)

public:

  virtual void getOptimalMapping( Graph<NodeAttribute,EdgeAttribute> * g1,
                                  Graph<NodeAttribute,EdgeAttribute> * g2,
                                  int * G1_to_G2, int * G2_to_G1 ) = 0;



  IPFPGraphEditDistanceMulti( EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction,
                              BipartiteGraphEditDistanceMulti<NodeAttribute,EdgeAttribute> * ed_multi,
                              int n_edit_paths
                            ):
    IPFPGraphEditDistance<NodeAttribute,EdgeAttribute>(costFunction, NULL),
    MultiGed<NodeAttribute,EdgeAttribute>(n_edit_paths),
    _ed_multi(ed_multi)
  {}


};

#endif

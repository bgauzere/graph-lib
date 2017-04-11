/**
 * @file IPFPZetaGraphEditDistance.h
 * @author Evariste <<evariste.daller@unicaen.fr>>
 * @version     0.0.1 - Mon Apr  10 2017
 *
 *
 * Description of the program objectives.
 * All necessary references.
 */

#ifndef __IPFPZETAGRAPHEDITDISTANCEMULTI_H__
#define __IPFPZETAGRAPHEDITDISTANCEMULTI_H__


#include "IPFPZetaGraphEditDistance.h"
#include "IPFPGraphEditDistanceMulti.h"



template<class NodeAttribute, class EdgeAttribute>
class IPFPZetaGraphEditDistanceMulti:
  public IPFPZetaGraphEditDistance<NodeAttribute, EdgeAttribute>,
  public IPFPGraphEditDistanceMulti<NodeAttribute, EdgeAttribute>
{

public: /* CONSTRUCTORS AND ACCESSORS */

  IPFPZetaGraphEditDistanceMulti( EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction,
                                  GraphEditDistance<NodeAttribute,EdgeAttribute> * ed_init,
                                  double zeta, int nep ):
    IPFPZetaGraphEditDistance<NodeAttribute,EdgeAttribute>(costFunction, ed_init, zeta),
    IPFPGraphEditDistanceMulti<NodeAttribute,EdgeAttribute>(costFunction, ed_init, nep)
  {}


  IPFPZetaGraphEditDistanceMulti( EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction,
                                  double zeta, int nep ):
    IPFPGraphEditDistance<NodeAttribute,EdgeAttribute>(costFunction, NULL),
    IPFPZetaGraphEditDistance<NodeAttribute,EdgeAttribute>(costFunction, zeta),
    IPFPGraphEditDistanceMulti<NodeAttribute,EdgeAttribute>(costFunction, NULL, nep)
  {}


};


#endif // __IPFPZETAGRAPHEDITDISTANCEMULTI_H__

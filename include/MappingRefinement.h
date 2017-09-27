/**
 * @file MappingRefinement.h
 * @author Evariste <<evariste.daller@unicarn.fr>>
 * @version Sep  25  2017
 */

#ifndef __MAPPINGREFINEMENT_H__
#define __MAPPINGREFINEMENT_H__

#include <graph.h>

/**
 * An algorithm refining a mapping from an initialization.
 *
 *   The given initialization is a mapping too, non optimal, 
 *   and a derived class of MappingRefinement will compute
 *   a better mapping from the init.
 */
template<class NodeAttribute, class EdgeAttribute>
  class MappingRefinement
{

 public:

  /**
   * Refine the mapping given in input as G1_to_G2 and G2_ti_G1
   * from g1 to g2 and returns a better mapping in these arrays
   */
  virtual void getBetterMapping( Graph<NodeAttribute, EdgeAttribute>* g1, Graph<NodeAttribute, EdgeAttribute>* g2,
					 int* G1_to_G2,  int* G2_to_G1, bool fromInit=false ) = 0;

  /**
   * Compute and return the cost of the given mapping from
   * g1 to g2.
   */
  virtual double mappingCost( Graph<NodeAttribute, EdgeAttribute>* g1, Graph<NodeAttribute, EdgeAttribute>* g2,
			      int* G1_to_G2,  int* G2_to_G1 ) = 0;
  
  
  virtual MappingRefinement * clone() const = 0;
  
};



#endif

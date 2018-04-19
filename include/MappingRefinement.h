/**
 * @file MappingRefinement.h
 * @author Evariste <<evariste.daller@unicarn.fr>>
 * @version Sep  25  2017
 */

#ifndef __MAPPINGREFINEMENT_H__
#define __MAPPINGREFINEMENT_H__

#include <list>
#include "graph.h"

/**
 * @brief An algorithm refining a mapping from an initialization.
 *
 *   The given initialization is a mapping too, non optimal,
 *   and a derived class of MappingRefinement will compute
 *   a better mapping from the init.
 */
template<class NodeAttribute, class EdgeAttribute>
  class MappingRefinement
{

 protected:

   bool _compute_equiv_mappings;
   std::list<unsigned int*> _equiv_G1toG2;
   std::list<unsigned int*> _equiv_G2toG1;

 public:

   bool compute_equiv_mappings(){ return _compute_equiv_mappings; }
   void compute_equiv_mappings(bool choice) { _compute_equiv_mappings = choice; }

   std::list<unsigned int*> & getEquivalentG1toG2() { return _equiv_G1toG2; }
   std::list<unsigned int*> & getEquivalentG2toG1() { return _equiv_G2toG1; }


  /**
   * @brief Refine the mapping given in input as G1_to_G2 and G2_ti_G1 from g1 to g2 and returns a better mapping in these arrays
   * @param G1_to_G2  The forward mapping (from g1 to g2) to refine
   * @param G2_to_G1  The corresponding reverse mapping (from g2 to g1). This parameter is useful with LSAPE mappings
   * @param fromInit  Allow to set up if the refined mapping should be computed from the given initialization or from a generated one
   */
  virtual void getBetterMapping( Graph<NodeAttribute, EdgeAttribute>* g1, Graph<NodeAttribute, EdgeAttribute>* g2,
			 	unsigned int* G1_to_G2,  unsigned int* G2_to_G1, bool fromInit=false ) = 0;

  /**
   * @brief Compute and return the cost of the given mapping from g1 to g2.
   */
  virtual double mappingCost( Graph<NodeAttribute, EdgeAttribute>* g1, Graph<NodeAttribute, EdgeAttribute>* g2,
			      unsigned int* G1_to_G2,  unsigned int* G2_to_G1 ) = 0;

  /**
   * @brief Clone the derivated object
   */
  virtual MappingRefinement * clone() const = 0;

  virtual ~MappingRefinement(){}

};



#endif

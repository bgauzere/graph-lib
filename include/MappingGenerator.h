/**
 * @file MappingGenerator.h
 * @author Evariste <<evariste.daller@unicaen.fr>>
 * @version  Sep  25  2017 
 */

#ifndef __MAPPING_GENERATOR_H
#define __MAPPING_GENERATOR_H

#include <list>


template<class NodeAttribute, class EdgeAttribute>
  class MappingGenerator {

 public:

  virtual std::list<int*> getMappings( Graph<NodeAttribute, EdgeAttribute>* g1, Graph<NodeAttribute, EdgeAttribute>* g2,
				       int k ) = 0;
  
};


#endif

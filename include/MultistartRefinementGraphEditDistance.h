/**
 * @file MultipleIPFPGraphEditDistance.h
 * @author Evariste <<evariste.daller@unicaen.fr>>
 * @version     Jun  9 2017
 *
 */

#ifndef __MULTISTARTREFINEMENTGED_H__
#define __MULTISTARTREFINEMENTGED_H__

#ifdef _OPENMP
  #include <omp.h>
#endif

#include <sys/time.h>
#include <list>
#include "GraphEditDistance.h"
#include "MultistartMappingRefinement.h"

#include "MultiGed.h"



/**
 * @brief Implements a MultistartMappingRefinement method adapted to the GED problem
 */
template<class NodeAttribute, class EdgeAttribute>
class MultistartRefinementGraphEditDistance:
  public virtual GraphEditDistance<NodeAttribute,EdgeAttribute>,
  public MultistartMappingRefinement<NodeAttribute, EdgeAttribute>
{

protected:

  MappingRefinement<NodeAttribute, EdgeAttribute> * method; //!< Storage of a predefined refinement method
  std::list<unsigned int*> refinedReverseMappings; //!< Storage of the reverse mappings needed for the (n+1)*(m+1) GED modelisation

  bool cleanMethod; //!< Delete the method in the destructor if true

public:


  /**
   * @brief Get the best possible mapping for the GED, with the given initialization and refinement methods
   */
  virtual void getOptimalMapping( Graph<NodeAttribute,EdgeAttribute> * g1,
                                  Graph<NodeAttribute,EdgeAttribute> * g2,
                                  unsigned int * G1_to_G2, unsigned int * G2_to_G1 );


  virtual void getBestMappingFromSet( MappingRefinement<NodeAttribute, EdgeAttribute> * algorithm,
                                      Graph<NodeAttribute,EdgeAttribute> * g1,
                                      Graph<NodeAttribute,EdgeAttribute> * g2,
                                      unsigned int * G1_to_G2, unsigned int * G2_to_G1,
                                      std::list<unsigned int*>& mappings );


  /**
   * First Constructor
   */
  MultistartRefinementGraphEditDistance( EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction,
                                 MappingGenerator<NodeAttribute,EdgeAttribute> * gen,
                                 int n_edit_paths,
                                 MappingRefinement<NodeAttribute, EdgeAttribute> * algorithm
                               ):
    GraphEditDistance<NodeAttribute,EdgeAttribute> (costFunction),
    MultistartMappingRefinement<NodeAttribute, EdgeAttribute> (gen, n_edit_paths),
    method(algorithm),
    cleanMethod(false)
  {}


  MultistartRefinementGraphEditDistance(
                        const MultistartRefinementGraphEditDistance<NodeAttribute, EdgeAttribute>& other
                        ):
    GraphEditDistance<NodeAttribute,EdgeAttribute> (other.cf),
    MultistartMappingRefinement<NodeAttribute, EdgeAttribute> (other.initGen->clone(), other.k),
    method(other.method->clone()),
    cleanMethod(true)
  {}


  ~MultistartRefinementGraphEditDistance(){
    if (cleanMethod){
      delete method;
      delete this->initGen;
    }
  }


  /**
   * @brief Returns refined mappings generated by the internal generator \ref initGen
   *
   *  The refinement method is \ref method
   */
  virtual const std::list<unsigned int*>&
  getBetterMappings( Graph<NodeAttribute,EdgeAttribute> * g1,
                     Graph<NodeAttribute,EdgeAttribute> * g2 );


  /**
   * @brief Returns refined mappings from the ones given in parameter, according to the method \method
   *
   *   Note that the returned mappings are the forward mappings of size (n+1). To get the
   *   (m+1) reverse mappings denoted by G2_to_G1, use the method \ref getReverseMappings.
   *
   * @param  g1          First graph
   * @param  g2          Second graph
   * @param  mapping     a list of arrays representing initial mappings
   * @note   Forward and reverse mappings are allocated on the heap and memory management is left to the user
   * @see getBetterMappings getReverseMappings
   */
  virtual const std::list<unsigned int*>&
  getBetterMappingsFromSet( Graph<NodeAttribute,EdgeAttribute> * g1,
                            Graph<NodeAttribute,EdgeAttribute> * g2,
                            std::list<unsigned int*>& mappings );


  /**
   * @brief Redefinition of \ref MultistartMappingRefinement::getBetterMappingsFromSet for the GED problem
   *
   *   Note that the returned mappings are the forward mappings of size (n+1). To get the
   *   (m+1) reverse mappings denoted by G2_to_G1, use the method \ref getReverseMappings.
   *
   * @param  algorithm   the refinement method
   * @param  g1          First graph
   * @param  g2          Second graph
   * @param  mapping     a list of arrays representing initial mappings
   * @note   Forward and reverse mappings are allocated on the heap and memory management is left to the user
   * @see getBetterMappings getReverseMappings
   */
  virtual const std::list<unsigned int*>&
  getBetterMappingsFromSet( MappingRefinement<NodeAttribute, EdgeAttribute> * algorithm,
                            Graph<NodeAttribute,EdgeAttribute> * g1,
                            Graph<NodeAttribute,EdgeAttribute> * g2,
                            std::list<unsigned int*>& mappings );

  /**
   * @brief Returns the last reverse mappings G2_to_G1 computed from \ref getBetterMappingsFromSet or \ref getBetterMappings
   */
  std::list<unsigned int*>& getReverseMappings(){ return refinedReverseMappings; }


  /**
   * Clone
   */
   virtual MultistartRefinementGraphEditDistance<NodeAttribute, EdgeAttribute>* clone() const {
     return new MultistartRefinementGraphEditDistance<NodeAttribute, EdgeAttribute>(*this);
   }

};

//---


template<class NodeAttribute, class EdgeAttribute>
void MultistartRefinementGraphEditDistance<NodeAttribute, EdgeAttribute>::
getOptimalMapping( Graph<NodeAttribute,EdgeAttribute> * g1,
                   Graph<NodeAttribute,EdgeAttribute> * g2,
                   unsigned int * G1_to_G2, unsigned int * G2_to_G1)
{
  this->getBestMapping(method, g1, g2, G1_to_G2, G2_to_G1);
}



template<class NodeAttribute, class EdgeAttribute>
void MultistartRefinementGraphEditDistance<NodeAttribute, EdgeAttribute>::
getBestMappingFromSet( MappingRefinement<NodeAttribute, EdgeAttribute> * algorithm,
                Graph<NodeAttribute,EdgeAttribute> * g1,
                Graph<NodeAttribute,EdgeAttribute> * g2,
                unsigned int * G1_to_G2, unsigned int * G2_to_G1,
                std::list<unsigned int*>& mappings )
{
  unsigned int n = g1->Size();
  unsigned int m = g2->Size();
  double cost = -1;

  // This will update refinedMappings and refinedReverseMappings
  getBetterMappingsFromSet(algorithm, g1, g2, mappings);

  // Look for the minimal cost mapping
  std::list<unsigned int*>::const_iterator itf_optim, itr_optim;
  for (std::list<unsigned int*>::const_iterator itf = this->refinedMappings.begin(), itr = this->refinedReverseMappings.begin();
       itf != this->refinedMappings.end() && itr != this->refinedReverseMappings.end();
       itf++, itr++ )
  {
    double current_cost =  algorithm->mappingCost(g1, g2, *itf, *itr);
    if (cost > current_cost || cost == -1){
         cost = current_cost;
         itf_optim = itf;
         itr_optim = itr;
    }
  }

  // Copy the optimal mapping
  for (unsigned int i=0; i<n; i++) G1_to_G2[i] = (*itf_optim)[i];
  for (unsigned int j=0; j<m; j++) G2_to_G1[j] = (*itr_optim)[j];

  // Remove the computed mappings (local to the object)
  for (std::list<unsigned int*>::const_iterator itf = this->refinedMappings.begin(), itr = this->refinedReverseMappings.begin();
       itf != this->refinedMappings.end() && itr != this->refinedReverseMappings.end();
       itf++, itr++ )
  {
    delete[] *itf;
    delete[] *itr;
  }
}




template<class NodeAttribute, class EdgeAttribute>
const std::list<unsigned int*>& MultistartRefinementGraphEditDistance<NodeAttribute, EdgeAttribute>::
getBetterMappings( Graph<NodeAttribute,EdgeAttribute> * g1,
                   Graph<NodeAttribute,EdgeAttribute> * g2 )
{
  std::list<unsigned int*> mappings = this->initGen->getMappings(g1, g2, this->k);
  const std::list<unsigned int*>& refined = this->getBetterMappingsFromSet(method, g1, g2, mappings);

  // Delete original (bipartite) mappings
  for (std::list<unsigned int*>::iterator it=mappings.begin();
       it != mappings.end();   it++)
    delete [] *it;

  return refined;
}




template<class NodeAttribute, class EdgeAttribute>
const std::list<unsigned int*>& MultistartRefinementGraphEditDistance<NodeAttribute, EdgeAttribute>::
getBetterMappingsFromSet( Graph<NodeAttribute,EdgeAttribute> * g1,
                          Graph<NodeAttribute,EdgeAttribute> * g2,
                          std::list<unsigned int*>& mappings )
{
  return getBetterMappingsFromSet( method, g1, g2, mappings);
}



template<class NodeAttribute, class EdgeAttribute>
const std::list<unsigned int*>& MultistartRefinementGraphEditDistance<NodeAttribute, EdgeAttribute>::
getBetterMappingsFromSet( MappingRefinement<NodeAttribute, EdgeAttribute> * algorithm,
                          Graph<NodeAttribute,EdgeAttribute> * g1,
                          Graph<NodeAttribute,EdgeAttribute> * g2,
                          std::list<unsigned int*>& mappings )
{
  unsigned int n = g1->Size();
  unsigned int m = g2->Size();

  this->refinedMappings.clear(); //G1_to_G2 in refinedMappings
  this->refinedReverseMappings.clear(); //G2_to_G1 in refinedReverseMappings

  // computation of LSAPE versions of the mappings
  std::list<unsigned int*> listG1toG2, listG2toG1;
  MultiGed<NodeAttribute, EdgeAttribute>::mappings_lsap2lsape(
    mappings, listG1toG2, listG2toG1, n, m
  );

  #pragma omp parallel
  {
    // Iterators (private) : mustn't be shared accross the threads
    typename std::list<unsigned int*>::iterator itf, itr;

   #pragma omp single nowait
   {
    for ( itf=listG1toG2.begin(), itr=listG2toG1.begin();
          itf != listG1toG2.end() && itr != listG2toG1.end();
          itr++, itf++ )
    {
      #pragma omp task untied
      {
        #ifdef _OPENMP
          MappingRefinement<NodeAttribute, EdgeAttribute> * local_method = algorithm->clone();
        #else
          MappingRefinement<NodeAttribute, EdgeAttribute> * local_method = algorithm;
        #endif

        local_method->getBetterMapping(g1, g2, *itf, *itr, true);
        #pragma omp critical
        {
          this->refinedMappings.insert(
            this->refinedMappings.end(),
            local_method->getEquivalentG1toG2().begin(),
            local_method->getEquivalentG1toG2().end()
          );
          this->refinedReverseMappings.insert(
            this->refinedReverseMappings.end(),
            local_method->getEquivalentG2toG1().begin(),
            local_method->getEquivalentG2toG1().end()
          );
        } // end omp critical

        #ifdef _OPENMP
          delete local_method;
        #endif
      } // end omp task
    } //end for
   } // end omp single
  } // end omp parallel

  return this->refinedMappings;
}
#endif // __MULTISTARTREFINEMENTGED_H__

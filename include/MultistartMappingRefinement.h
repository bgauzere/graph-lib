/**
 * @file MultipleIPFPGraphEditDistance.h
 * @author Evariste <<evariste.daller@unicaen.fr>>
 * @version     Jun  9 2017
 *
 */

#ifndef __MULTISTARTMAPPINGREFINEMENT_H__
#define __MULTISTARTMAPPINGREFINEMENT_H__

#ifdef _OPENMP
  #include <omp.h>
#endif

#include <sys/time.h>
#include <list>
#include "MappingRefinement.h"
#include "MappingGenerator.h"


template<class NodeAttribute, class EdgeAttribute>
class MultistartMappingRefinement
{

protected:

  MappingGenerator<NodeAttribute, EdgeAttribute> * initGen; //! Generator of initializations 
  int k;

public:


  virtual void getBestMapping( MappingRefinement<NodeAttribute, EdgeAttribute> * algorithm,
                               Graph<NodeAttribute,EdgeAttribute> * g1,
                               Graph<NodeAttribute,EdgeAttribute> * g2,
                               int * G1_to_G2, int * G2_to_G1 );

  virtual void getBestMappingFromSet( MappingRefinement<NodeAttribute, EdgeAttribute> * algorithm,
                                      Graph<NodeAttribute,EdgeAttribute> * g1,
                                      Graph<NodeAttribute,EdgeAttribute> * g2,
                                      int * G1_to_G2, int * G2_to_G1,
                                      std::list<int*>& mappings );


  MultistartMappingRefinement( MappingGenerator<NodeAttribute, EdgeAttribute> * gen,
                               int nSol
                             ):
    initGen(gen),
    k(nSol)
  {}


};

//---


template<class NodeAttribute, class EdgeAttribute>
void MultistartMappingRefinement<NodeAttribute, EdgeAttribute>::
getBestMapping( MappingRefinement<NodeAttribute, EdgeAttribute> * algorithm,
                Graph<NodeAttribute,EdgeAttribute> * g1,
                Graph<NodeAttribute,EdgeAttribute> * g2,
                int * G1_to_G2, int * G2_to_G1 )
{
  //Compute Mapping init
  struct timeval  tv1, tv2;
  gettimeofday(&tv1, NULL);

  std::list<int*> mappings = initGen->getMappings(g1, g2, k);
  gettimeofday(&tv2, NULL);

  this->getBestMappingFromSet(algorithm, g1, g2, G1_to_G2, G2_to_G1, mappings);
}



template<class NodeAttribute, class EdgeAttribute>
void MultistartMappingRefinement<NodeAttribute, EdgeAttribute>::
getBestMappingFromSet( MappingRefinement<NodeAttribute, EdgeAttribute> * algorithm,
                Graph<NodeAttribute,EdgeAttribute> * g1,
                Graph<NodeAttribute,EdgeAttribute> * g2,
                int * G1_to_G2, int * G2_to_G1,
                std::list<int*>& mappings )
{
  struct timeval  tv1, tv2;
  int n = g1->Size();
  int m = g2->Size();

  typename std::list<int*>::const_iterator it;
  double cost = -1;
  double ncost;


  // Multithread
  #ifdef _OPENMP
    gettimeofday(&tv1, NULL);
    int** arrayMappings = new int*[mappings.size()];
    int* arrayCosts = new int[mappings.size()];
    int* arrayLocal_G1_to_G2 = new int[n * mappings.size()];

    int i=0; for (it=mappings.begin(); it!=mappings.end(); it++){
      arrayMappings[i] = *it;
      i++;
    }

    //omp_set_dynamic(0);
    omp_set_num_threads(4);
    #pragma omp parallel for schedule(dynamic) //private(tid, i, j, ncost, ipfpGed )
    for (int tid=0; tid<mappings.size(); tid++){
      int* lsapMapping = arrayMappings[tid];
      int* local_G1_to_G2 = &(arrayLocal_G1_to_G2[tid*n]);

  // Sequential
  #else
    int* local_G1_to_G2 = new int[n];

    double t_acc = 0; // accumulated time
    for (it=mappings.begin(); it!=mappings.end(); it++){
      gettimeofday(&tv1, NULL);
      int* lsapMapping = *it;
  #endif

    // Copy the mapping into the local array
    for (int i=0; i<n; i++)
      local_G1_to_G2[i] = lsapMapping[i];
    
    MappingRefinement<NodeAttribute, EdgeAttribute> * local_method;
    
    #ifdef _OPENMP
      local_method = algorithm->clone();
    #else
      local_method = algorithm;
    #endif
    
    local_method->getBetterMapping(g1, g2, local_G1_to_G2, NULL, true);
    ncost = local_method->mappingCost(g1, g2, local_G1_to_G2, NULL);


    // Multithread
    #ifdef _OPENMP
      // save the approx cost
      arrayCosts[tid] = ncost;
     // _distances_[tid] = ncost;

    // Sequential
    #else
      // if ncost is better : save the mapping and the cost
      if (cost > ncost || cost == -1){
        cost = ncost;
        for (int i=0; i<n; i++) G1_to_G2[i] = local_G1_to_G2[i];
      }
      gettimeofday(&tv2, NULL);
      t_acc += ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));

    #endif

  } //end for

  // Multithread : Reduction
  #ifdef _OPENMP
    gettimeofday(&tv2, NULL);
    
    gettimeofday(&tv1, NULL);

    int i_optim;
    for (int i=0; i<mappings.size(); i++){
      if (cost > arrayCosts[i] || cost == -1){
         cost = arrayCosts[i];
         i_optim = i;
      }
    }
    for (int i=0; i<n; i++) G1_to_G2[i] = arrayLocal_G1_to_G2[i_optim*n + i];

    // To match the output format size in XPs
    //for (int i=mappings.size(); i<k; i++) _distances_[i] = 9999;

    gettimeofday(&tv2, NULL);
    //_xp_out_ <<  ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec)) << ", ";
    //_xp_out_ << ((float)t) / CLOCKS_PER_SEC << ", ";

    delete[] arrayLocal_G1_to_G2;
    delete[] arrayCosts;
    delete[] arrayMappings;

  // Sequential : deletes
  #else

    delete [] local_G1_to_G2;

  #endif


}


#endif // __MULTIPLEIPFPGRAPHEDITDISTANCE_H__

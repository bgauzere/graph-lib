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


template<class NodeAttribute, class EdgeAttribute>
class MultistartRefinementGraphEditDistance:
  public virtual GraphEditDistance<NodeAttribute,EdgeAttribute>,
  public MultistartMappingRefinement<NodeAttribute, EdgeAttribute>
{

protected:

  MappingRefinement<NodeAttribute, EdgeAttribute> * method;

public:


  virtual void getOptimalMapping( Graph<NodeAttribute,EdgeAttribute> * g1,
                                  Graph<NodeAttribute,EdgeAttribute> * g2,
                                  int * G1_to_G2, int * G2_to_G1 );


  virtual void getBestMappingFromSet( MappingRefinement<NodeAttribute, EdgeAttribute> * algorithm,
                                      Graph<NodeAttribute,EdgeAttribute> * g1,
                                      Graph<NodeAttribute,EdgeAttribute> * g2,
                                      int * G1_to_G2, int * G2_to_G1,
                                      std::list<int*>& mappings );


  MultistartRefinementGraphEditDistance( EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction,
                                 MappingGenerator<NodeAttribute,EdgeAttribute> * gen,
                                 int n_edit_paths,
                                 MappingRefinement<NodeAttribute, EdgeAttribute> * algorithm
                               ):
    GraphEditDistance<NodeAttribute,EdgeAttribute> (costFunction),
    MultistartMappingRefinement<NodeAttribute, EdgeAttribute> (gen, n_edit_paths),
    method(algorithm)
  {}


};

//---


template<class NodeAttribute, class EdgeAttribute>
void MultistartRefinementGraphEditDistance<NodeAttribute, EdgeAttribute>::
getOptimalMapping( Graph<NodeAttribute,EdgeAttribute> * g1,
                   Graph<NodeAttribute,EdgeAttribute> * g2,
                   int * G1_to_G2, int * G2_to_G1)
{
  this->getBestMapping(method, g1, g2, G1_to_G2, G2_to_G1);
}



template<class NodeAttribute, class EdgeAttribute>
void MultistartRefinementGraphEditDistance<NodeAttribute, EdgeAttribute>::
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
    int* arrayLocal_G1_to_G2 = new int[(n+1) * mappings.size()];
    int* arrayLocal_G2_to_G1 = new int[(m+1) * mappings.size()];

    int i=0; for (it=mappings.begin(); it!=mappings.end(); it++){
      arrayMappings[i] = *it;
      i++;
    }

    //omp_set_dynamic(0);
    omp_set_num_threads(4);
    #pragma omp parallel for schedule(dynamic) //private(tid, i, j, ncost, ipfpGed )
    for (int tid=0; tid<mappings.size(); tid++){
      int* lsapMapping = arrayMappings[tid];
      int* local_G1_to_G2 = &(arrayLocal_G1_to_G2[tid*(n+1)]);
      int* local_G2_to_G1 = &(arrayLocal_G2_to_G1[tid*(m+1)]);

  // Sequential
  #else
    int* local_G1_to_G2 = new int[n+1];
    int* local_G2_to_G1 = new int[m+1];

    double t_acc = 0; // accumulated time
    for (it=mappings.begin(); it!=mappings.end(); it++){
      gettimeofday(&tv1, NULL);
      int* lsapMapping = *it;
  #endif

    // computation of G1_to_G2 and G2_to_G1
    for (int j=0; j<m; j++){ // connect all to epsilon by default
      local_G2_to_G1[j] = n;
    }

    for (int i=0; i<n; i++){
      if (lsapMapping[i] >= m)
        local_G1_to_G2[i] = m; // i -> epsilon
      else{
        local_G1_to_G2[i] = lsapMapping[i];
        local_G2_to_G1[lsapMapping[i]] = i;
      }
    }

    for (int j=0; j<m; j++){
      if (lsapMapping[n+j] < m){
        local_G2_to_G1[j] = n; // epsilon -> j
      }
    }
    
    MappingRefinement<NodeAttribute, EdgeAttribute> * local_method;
    
    #ifdef _OPENMP
      local_method = algorithm->clone();
    #else
      local_method = algorithm;
    #endif
    
    local_method->getBetterMapping(g1, g2, local_G1_to_G2, local_G2_to_G1);
    ncost = local_method->mappingCost(g1, g2, local_G1_to_G2, local_G2_to_G1);


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
        for (int j=0; j<m; j++) G2_to_G1[j] = local_G2_to_G1[j];
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
    for (int i=0; i<n; i++) G1_to_G2[i] = arrayLocal_G1_to_G2[i_optim*(n+1)+i];
    for (int j=0; j<m; j++) G2_to_G1[j] = arrayLocal_G2_to_G1[i_optim*(m+1)+j];

    // To match the output format size in XPs
    //for (int i=mappings.size(); i<k; i++) _distances_[i] = 9999;

    gettimeofday(&tv2, NULL);
    //_xp_out_ <<  ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec)) << ", ";
    //_xp_out_ << ((float)t) / CLOCKS_PER_SEC << ", ";

    delete[] arrayLocal_G1_to_G2;
    delete[] arrayLocal_G2_to_G1;
    delete[] arrayCosts;
    delete[] arrayMappings;

  // Sequential : deletes
  #else

    delete [] local_G1_to_G2;
    delete [] local_G2_to_G1;

  #endif

}

#endif // __MULTISTARTREFINEMENTGED_H__

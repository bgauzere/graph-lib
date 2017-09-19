/**
 * @file MultipleIPFPGraphEditDistance.h
 * @author Evariste <<evariste.daller@unicaen.fr>>
 * @version     Jun  9 2017
 *
 */

#ifndef __MULTIPLEIPFPGRAPHEDITDISTANCE_H__
#define __MULTIPLEIPFPGRAPHEDITDISTANCE_H__

#ifdef _OPENMP
  #include <omp.h>
#endif

#include <sys/time.h>
#include <list>
#include "GraphEditDistance.h"
#include "IPFPGraphEditDistance.h"
#include "BipartiteGraphEditDistanceMulti.h"


template<class NodeAttribute, class EdgeAttribute>
class MultipleIPFPGraphEditDistance:
  public GraphEditDistance<NodeAttribute,EdgeAttribute>
{

protected:

  BipartiteGraphEditDistanceMulti<NodeAttribute,EdgeAttribute> * ed_multi;  //! The multiple solution based GED used (ie. Bunke or Random Walks)
  int nep;

public:


  virtual void getOptimalMapping( Graph<NodeAttribute,EdgeAttribute> * g1,
                                  Graph<NodeAttribute,EdgeAttribute> * g2,
                                  int * G1_to_G2, int * G2_to_G1 );



  MultipleIPFPGraphEditDistance( EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction,
                                 BipartiteGraphEditDistanceMulti<NodeAttribute,EdgeAttribute> * ed_multi,
                                 int n_edit_paths
                               ):
    GraphEditDistance<NodeAttribute,EdgeAttribute> (costFunction),
    ed_multi(ed_multi),
    nep(n_edit_paths)
  {}


};

//---


template<class NodeAttribute, class EdgeAttribute>
void MultipleIPFPGraphEditDistance<NodeAttribute, EdgeAttribute>::
getOptimalMapping( Graph<NodeAttribute,EdgeAttribute> * g1,
                   Graph<NodeAttribute,EdgeAttribute> * g2,
                   int * G1_to_G2, int * G2_to_G1)
{
  //Compute Mapping init
  struct timeval  tv1, tv2;
  gettimeofday(&tv1, NULL);

  std::list<int*> mappings = ed_multi->getOptimalMappings(g1, g2, nep);
  gettimeofday(&tv2, NULL);


  int n = g1->Size();
  int m = g2->Size();

  typename std::list<int*>::const_iterator it;
  double ged = -1;
  double nged;


  // Multithread
  #ifdef _OPENMP
    gettimeofday(&tv1, NULL);
    int** arrayMappings = new int*[mappings.size()];
    int* arrayGed = new int[mappings.size()];
    int* arrayLocal_G1_to_G2 = new int[(n+1) * mappings.size()];
    int* arrayLocal_G2_to_G1 = new int[(m+1) * mappings.size()];

    int i=0; for (it=mappings.begin(); it!=mappings.end(); it++){
      arrayMappings[i] = *it;
      i++;
    }

    //omp_set_dynamic(0);
    omp_set_num_threads(4);
    #pragma omp parallel for schedule(dynamic) //private(tid, i, j, nged, ipfpGed )
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

    IPFPGraphEditDistance<NodeAttribute,EdgeAttribute> ipfpGed(this->cf);
    ipfpGed.getOptimalMappingFromInit(g1, g2, local_G1_to_G2, local_G2_to_G1);
    //IPFPQAP<NodeAttribute, EdgeAttribute> ipfpGed(this->cf);
    //ipfpGed.getOptimalMapping(g1, g2, local_G1_to_G2, local_G2_to_G1);

    nged = this->GedFromMapping(g1, g2, local_G1_to_G2, n, local_G2_to_G1, m);

    // Multithread
    #ifdef _OPENMP
      // save the approx ged
      arrayGed[tid] = nged;
      //_distances_[tid] = nged;

    // Sequential
    #else
      // if nged is better : save the mapping and the ged
      if (ged > nged || ged == -1){
        ged = nged;
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
    //_xp_out_ <<  ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec)) << ", ";
    //_xp_out_ << ", " << ((float)t) / CLOCKS_PER_SEC;
    gettimeofday(&tv1, NULL);

    int i_optim;
    for (int i=0; i<mappings.size(); i++){
      if (ged > arrayGed[i] || ged == -1){
         ged = arrayGed[i];
         i_optim = i;
      }
    }
    for (int i=0; i<n; i++) G1_to_G2[i] = arrayLocal_G1_to_G2[i_optim*(n+1)+i];
    for (int j=0; j<m; j++) G2_to_G1[j] = arrayLocal_G2_to_G1[i_optim*(m+1)+j];

    // To match the output format size in XPs
    //for (int i=mappings.size(); i<nep; i++) _distances_[i] = 9999;

    gettimeofday(&tv2, NULL);
    //_xp_out_ <<  ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec)) << ", ";
    //_xp_out_ << ((float)t) / CLOCKS_PER_SEC << ", ";

    delete[] arrayLocal_G1_to_G2;
    delete[] arrayLocal_G2_to_G1;
    delete[] arrayGed;
    delete[] arrayMappings;

  // Sequential : deletes
  #else

    delete [] local_G1_to_G2;
    delete [] local_G2_to_G1;
  #endif


}


#endif // __MULTIPLEIPFPGRAPHEDITDISTANCE_H__

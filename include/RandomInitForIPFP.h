#ifndef __RANDFORIPFP_H__
#define __RANDFORIPFP_H__


#include "BipartiteGraphEditDistanceMulti.h"

template<class NodeAttribute, class EdgeAttribute>
class RandomInitForIPFP :
  public BipartiteGraphEditDistanceMulti<NodeAttribute, EdgeAttribute>
{

public:

  RandomInitForIPFP( EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction,
                     int _k ):
    BipartiteGraphEditDistance<NodeAttribute,EdgeAttribute>(costFunction),
    BipartiteGraphEditDistanceMulti<NodeAttribute,EdgeAttribute>(costFunction, _k)
  {};

public:

  virtual std::list<int*> getOptimalMappings( Graph<NodeAttribute,EdgeAttribute> * g1,
                                              Graph<NodeAttribute,EdgeAttribute> * g2,
                                              int k = -1 );

  void getOptimalMapping( Graph<NodeAttribute,EdgeAttribute> * g1,
                          Graph<NodeAttribute,EdgeAttribute> * g2,
                          int * G1_to_G2, int * G2_to_G1 );

};


template<class NodeAttribute, class EdgeAttribute>
std::list<int*> RandomInitForIPFP<NodeAttribute, EdgeAttribute>::
getOptimalMappings( Graph<NodeAttribute,EdgeAttribute> * g1,
                    Graph<NodeAttribute,EdgeAttribute> * g2,
                    int k )
{
  if (k < 0) k = 100;
  
  int n = g1->Size();
  int m = g2->Size();
  
  std::list<int*> mappings;
  
  for (int i=0; i<k; i++){
    int* _map = new int[n+m];
    for (int a=0; a<n+m; a++) _map[a] = a;
    
    std::random_shuffle(&_map[0], &_map[n+m]);
    mappings.push_back(_map);
  }
  
  return mappings;
}


template<class NodeAttribute, class EdgeAttribute>
void RandomInitForIPFP<NodeAttribute, EdgeAttribute>::
getOptimalMapping (Graph<NodeAttribute,EdgeAttribute> * g1,
                   Graph<NodeAttribute,EdgeAttribute> * g2,
                   int * G1_to_G2, int * G2_to_G1 )
{
  std::list<int*> mappings = getOptimalMappings(g1, g2, this->_nep);
  
  typename std::list<int*>::const_iterator it;


  // Get the min of ged;
  this->_ged = -1;
  double nged;
  int n=g1->Size();
  int m=g2->Size();

  int* local_G1_to_G2 = new int[n+1];
  int* local_G2_to_G1 = new int[m+1];

  for (it=mappings.begin(); it!=mappings.end(); it++){
    int* lsapMapping = *it;

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

    nged = this->GedFromMapping(g1, g2, local_G1_to_G2,n, local_G2_to_G1,m);


    // if nged is better : save the mapping and the ged
    if (this->_ged > nged || this->_ged == -1){
      this->_ged = nged;
      for (int i=0; i<n; i++) G1_to_G2[i] = local_G1_to_G2[i];
      for (int j=0; j<m; j++) G2_to_G1[j] = local_G2_to_G1[j];
    }
  }

  delete [] local_G1_to_G2;
  delete [] local_G2_to_G1;


  typename std::list<int*>::iterator it_del;
  for (it_del=mappings.begin(); it_del!=mappings.end(); it_del++)
    delete[] *it_del;
}

#endif

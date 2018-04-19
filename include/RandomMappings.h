#ifndef __RANDOMMAPPINGS_H__
#define __RANDOMMAPPINGS_H__
#include <random>
#include "MappingGenerator.h"

template<class NodeAttribute, class EdgeAttribute>
class RandomMappings :
  public MappingGenerator<NodeAttribute, EdgeAttribute>
{
protected:
  int _k;
  std::default_random_engine randGen;
public:

  RandomMappings(int k, unsigned int seed = 123):_k(k)
  {
    randGen.seed(seed);
  }
  
  virtual ~RandomMappings(){}

public:

  virtual std::list<unsigned int*> getMappings( Graph<NodeAttribute,EdgeAttribute> * g1,
						Graph<NodeAttribute,EdgeAttribute> * g2);


  virtual RandomMappings<NodeAttribute, EdgeAttribute> * clone() const {
   return new RandomMappings<NodeAttribute, EdgeAttribute>(*this);
  }
};



template<class NodeAttribute, class EdgeAttribute>
class RandomMappingsGED :
  public RandomMappings<NodeAttribute, EdgeAttribute>
{

public:

  RandomMappingsGED(int k, unsigned int seed = 123) :
    RandomMappings<NodeAttribute, EdgeAttribute>(k, seed)
  {}
  
  ~RandomMappingsGED(){}

public:

  virtual std::list<unsigned int*> getMappings( Graph<NodeAttribute,EdgeAttribute> * g1,
						Graph<NodeAttribute,EdgeAttribute> * g2);

};





template<class NodeAttribute, class EdgeAttribute>
std::list<unsigned int*> RandomMappings<NodeAttribute, EdgeAttribute>::
getMappings( Graph<NodeAttribute,EdgeAttribute> * g1,
	     Graph<NodeAttribute,EdgeAttribute> * g2)
{
  if (this->_k < 0) this->_k = 100;

  int n = g1->Size();
  int m = g2->Size();

  std::list<unsigned int*> mappings;

  int mx = std::max(n,m);

  for (int i=0; i<this->_k; i++){
    unsigned int* _map = new unsigned int[mx];
    for (int a=0; a<m; a++)  _map[a] = a;
    for (int a=m; a<mx; a++) _map[a] = -1; // if m<max(n,m) : complete with -1

    std::shuffle(&_map[0], &_map[mx], randGen);


    if (mx != n){
      unsigned int* mapping = new unsigned int[n];
      memcpy(mapping, _map, n*sizeof(int));
      mappings.push_back(mapping);
      delete [] _map;
    }
    else
      mappings.push_back(_map);
  }

  return mappings;
}




template<class NodeAttribute, class EdgeAttribute>
std::list<unsigned int*> RandomMappingsGED<NodeAttribute, EdgeAttribute>::
getMappings( Graph<NodeAttribute,EdgeAttribute> * g1,
	     Graph<NodeAttribute,EdgeAttribute> * g2)
{
  if (this->_k < 0) this->_k = 100;
  
  int n = g1->Size();
  int m = g2->Size();
  
  std::list<unsigned int*> mappings;
  
  for (int i=0; i<this->_k; i++){
    unsigned int* _map = new unsigned int[n+m];
    for (int a=0; a<n+m; a++) _map[a] = a;
    
    std::shuffle(&_map[0], &_map[n+m], this->randGen);
    mappings.push_back(_map);
  }
  
  return mappings;
}

#endif

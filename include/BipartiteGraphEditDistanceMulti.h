/**
 * @file BipartiteGraphEditDistanceMulti.h
 * @author Évariste <<evariste.daller@unicaen.fr>>
 * @version     0.0.1 - Mon Mar  27 2017
 *
 */

#ifndef __BIPARTITEGRAPHEDITDISTANCEMULTI_H__
#define __BIPARTITEGRAPHEDITDISTANCEMULTI_H__


#include "lsape.h"
#include "BipartiteGraphEditDistance.h"


/**
 * @brief Bipartite version of MultiGed based on the Bunke method for linear cost matrix generation
 */
template<class NodeAttribute, class EdgeAttribute>
class BipartiteGraphEditDistanceMulti:
  public BipartiteGraphEditDistance<NodeAttribute, EdgeAttribute>{
private:
  int _k;

public:
  
  BipartiteGraphEditDistanceMulti(EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction,
				  int _k ):
    BipartiteGraphEditDistance<NodeAttribute,EdgeAttribute>(costFunction),_k(_k)
  {};
  
  
  BipartiteGraphEditDistanceMulti(const BipartiteGraphEditDistanceMulti<NodeAttribute,EdgeAttribute> & other):
    BipartiteGraphEditDistance<NodeAttribute,EdgeAttribute>(other.cf),_k(other._k)   
  {
    this->C = NULL;
  }


public:

  /**
   * @brief compute an optimal mapping between `g1` and `g2`
   *        from k different optimal mappings by minimizing the ged optained
   * @note The GED is computed and set in `ged`
   */
  virtual void getOptimalMapping( Graph<NodeAttribute,EdgeAttribute> * g1,
                                  Graph<NodeAttribute,EdgeAttribute> * g2,
                                  int * G1_to_G2, int * G2_to_G1 );
  virtual ~BipartiteGraphEditDistanceMulti(){
    if (this->C != NULL) {
      delete [] this->C;
      this->C = NULL;
    }
  }
};

/* Optimal mapping is computed according to best obtained GED among all mappings*/
template<class NodeAttribute, class EdgeAttribute>
void BipartiteGraphEditDistanceMulti<NodeAttribute, EdgeAttribute>::
getOptimalMapping (Graph<NodeAttribute,EdgeAttribute> * g1,
                   Graph<NodeAttribute,EdgeAttribute> * g2,
                   int * G1_to_G2, int * G2_to_G1 )
{
  int n=g1->Size();
  int m=g2->Size();
  // Compute C
  delete [] this->C;
  this->computeCostMatrix(g1,g2);
  //Compute optimal assignement
  double *u = new double[n+1];
  double *v = new double[m+1];
  std::list<unsigned int*> solutions;

  __attribute__((unused)) double min_cost = lsape::lsapeSolutions<double>(this->C,n+1,m+1,this->_k,
						  solutions, this->my_solver);
  delete [] u;
  delete [] v;
  typename std::list<unsigned int*>::const_iterator it;
  
  // Get the min of ged;
  int _ged = -1;
  double nged;

  unsigned int* local_G1_to_G2 = new unsigned int[n+1];
  unsigned int* local_G2_to_G1 = new unsigned int[m+1];

  for (it=solutions.begin(); it!=solutions.end(); it++){
    unsigned int* lsapMapping = *it;
    
    // computation of G1_to_G2 and G2_to_G1
    for (int j=0; j<m; j++){ // connect all to epsilon by default
      local_G2_to_G1[j] = n;
    }

    for (int i=0; i<n; i++){
      if (lsapMapping[i] >= (unsigned int)(m))
        local_G1_to_G2[i] = m; // i -> epsilon
      else{
        local_G1_to_G2[i] = lsapMapping[i];
        local_G2_to_G1[lsapMapping[i]] = i;
      }
    }

    for (int j=0; j<m; j++){
      if (lsapMapping[n+j] <  (unsigned int)(m)){
        local_G2_to_G1[j] = n; // epsilon -> j
      }
    }

    nged = this->GedFromMapping(g1, g2, local_G1_to_G2,n, local_G2_to_G1,m);

    // if nged is better : save the mapping and the ged
    if (_ged > nged || _ged == -1){
      _ged = nged;
      for (int i=0; i<n; i++) G1_to_G2[i] = local_G1_to_G2[i];
      for (int j=0; j<m; j++) G2_to_G1[j] = local_G2_to_G1[j];
    }
  }

  delete [] local_G1_to_G2;
  delete [] local_G2_to_G1;

  typename std::list<unsigned int*>::iterator it_del;
  for (it_del=solutions.begin(); it_del!=solutions.end(); it_del++)
    delete[] *it_del;
}

#endif

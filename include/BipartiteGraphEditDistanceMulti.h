/**
 * @file BipartiteGraphEditDistanceMulti.h
 * @author Évariste <<evariste.daller@unicaen.fr>>
 * @version     0.0.1 - Mon Mar  27 2017
 */

#ifndef __BIPARTITEGRAPHEDITDISTANCEMULTI_H__
#define __BIPARTITEGRAPHEDITDISTANCEMULTI_H__

#include "BipartiteGraphEditDistance.h"
#include "AllPerfectMatchings.h"

template<class NodeAttribute, class EdgeAttribute>
class BipartiteGraphEditDistanceMulti :
      public BipartiteGraphEditDistance<NodeAttribute, EdgeAttribute>
{
private:
  int k;

public:

  BipartiteGraphEditDistanceMulti( EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction,
                                   int _k ):
    BipartiteGraphEditDistance<NodeAttribute,EdgeAttribute>(costFunction),
    k(_k)
  {};

  void setK(int newk) { k = newk; } //!< Pre-set k before runing operator()(Graph*, Graph*)
  int getK() { return k; }

public:

  /**
   * @brief Allocate and retruns $k$ optimal mappings between <code>g1</code> and <code>g2</code>
   * @param k  The number of mappings to compute, -1 to get all perfect matchings
   * @return  A list of mappings given as arrays of int. For each mapping M, <code>M[i]</code> is the mapping, in g2, of node i in g1
   * @note  Each array is allocated here and have to be deleted manually
   */
  virtual std::list<int*> getKOptimalMappings(Graph<NodeAttribute,EdgeAttribute> * g1,
                                              Graph<NodeAttribute,EdgeAttribute> * g2,
                                              const int& k=-1);

  /**
   * @brief Compute the Graph Edit Distance between <code>g1</code> and <code>g2</code> considering $k$ edit paths
   * @param k  The number of edit paths to compute
   */
  virtual double operator() (Graph<NodeAttribute,EdgeAttribute> * g1,
   			                     Graph<NodeAttribute,EdgeAttribute> * g2,
                             const int& k=-1);

  /**
   * @brief Compute the GED between <code>g1</code> and <code>g2</code> as the minimum GED found trough all edit paths
   * @return  calls to operator() (g1, g2, k=1)
   */
  virtual double operator() (Graph<NodeAttribute,EdgeAttribute> * g1,
                             Graph<NodeAttribute,EdgeAttribute> * g2);
};


//---


template<class NodeAttribute, class EdgeAttribute>
std::list<int*> BipartiteGraphEditDistanceMulti<NodeAttribute, EdgeAttribute>::
getKOptimalMappings(Graph<NodeAttribute,EdgeAttribute> * g1,
                    Graph<NodeAttribute,EdgeAttribute> * g2,
                    const int& k)
{
  int n=g1->Size();
  int m=g2->Size();

  // Compute C
  delete [] this->C;
  this->computeCostMatrix(g1, g2);

  // Compute an optimal assignement
  double *u = new double[n+1];
  double *v = new double[m+1];
  int* G1_to_G2 = new int[n+1];
  int* G2_to_G1 = new int[m+1];
  hungarianLSAPE(this->C,n+1,m+1, G1_to_G2, G2_to_G1, u,v,false);

  // Compute the k optimal mappings
  cDigraph<int> edg = equalityDigraph<double,int> (this->C, n, m, G1_to_G2, u, v);
  AllPerfectMatchings<int> apm(edg);
  apm.enumPerfectMatchings(edg,k);
  std::list<int*> mappings = apm.getPerfectMatchings();

  // Add the first one to the list
  mappings.push_front(G1_to_G2);

  delete[] G2_to_G1;
  delete [] u;
  delete [] v;

  return mappings;
}



template<class NodeAttribute, class EdgeAttribute>
double BipartiteGraphEditDistanceMulti<NodeAttribute, EdgeAttribute>::
operator() (Graph<NodeAttribute,EdgeAttribute> * g1,
            Graph<NodeAttribute,EdgeAttribute> * g2,
            const int& k)
{
  int n=g1->Size();
  int m=g2->Size();

  std::list<int*> mappings = getKOptimalMappings(g1, g2, k);

  // Get the min of ged;
  double ged = -1;  double nged;
  int* G2_to_G1 = new int[m+1];

  typename std::list<int*>::const_iterator it;
  for (it=mappings.begin(); it!=mappings.end(); it++){
    int* G1_to_G2 = *it;

    // computation of G2_to_G1
    for (int j=0; j<m; j++) G2_to_G1[j] = 0; // deconnect all
    for (int i=0; i<n; i++) G2_to_G1[G1_to_G2[i]] = i; // construct G2_to_G1
    nged = this->GedFromMapping(g1, g2, G1_to_G2,n, G2_to_G1,m);

    if (ged > nged || ged == -1)
      ged = nged;
  }

  for (it=mappings.begin(); it!=mappings.end(); it++)
    delete[] *it;

  delete [] G2_to_G1;
  return ged;
}


template<class NodeAttribute, class EdgeAttribute>
double BipartiteGraphEditDistanceMulti<NodeAttribute, EdgeAttribute>::
operator() (Graph<NodeAttribute,EdgeAttribute> * g1,
            Graph<NodeAttribute,EdgeAttribute> * g2)
{
  return operator() (g1,g2,k);
}

#endif

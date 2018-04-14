/**
 * @file MultiGed.h
 * @author Évariste <<evariste.daller@unicaen.fr>>
 * @version     0.0.1 -  Apr 6 2017
 * Given a cost matrix, compute the best mapping according to a given GraphEditDistance
 *
 */

#ifndef __MULTIGED_H__
#define __MULTIGED_H__

#include "GraphEditDistance.h"
#include <lsap.h>
#include <lsape.h>
#include "MappingGenerator.h"
#include "gl_utils.h"

/**
 * @brief Base class for graph edit distance methods that consider multiple linear approximations via a modified Uno algorithm
 *
 *   This class implements some methods to compute several assignments through the LSAPE library.
 *   It may be considered as a MappingGenerator using Uno's algorithm to compute bipartite optimal mappings
 *   from an arbitrary one, solution of the hungarian algorithm.
 *
 * @see RandomMappings
 * @see GreedyGraphEditDistance
 */
template<class NodeAttribute, class EdgeAttribute>
  class MultiGed{

protected: /* MEMBERS */

  int _nep; //!< number of edit paths to compute GED
  double _ged;


public: /* STATIC FUNCTIONS */

  /**
   * @brief Compute the cost matrix corresponding to the lsap given the lsape matrix C
   * @return  A heap allocated double matrix of size (n+m)*(n+m)
   */
  static double* costMatrixLSAP(double* C, int n, int m);


  /**
   * @brief Allocate and retruns `k` optimal mappings between `g1` and `g2`
   * @param k  The number of mappings to compute, -1 to get all perfect matchings
   * @param C  The cost matrix of size \f$(n+1)\times (m+1)\f$
   * @param n,m  size of the graphs
   * @return  A list of mappings given as arrays of int. For each mapping M, `M[i]` is the mapping, in g2, of node i in g1
   * @note  Each array is allocated here and have to be deleted manually
   */
  static std::list<unsigned int*> getKOptimalMappings(double* C,
						      int n, int m,
						      const int& k);
public: /* CONSTRUCTORS AND ACCESSORS */

  MultiGed(int _k):
    _nep(_k), _ged(-1)
  {}


  virtual ~MultiGed(){}


  virtual void setK(int newk) { _nep = newk; } //!< Set the number of assignments to be calculated to newk
  virtual int getK() { return _nep; }  //!< Returns the number of assignments

  virtual double getGED() {return _ged; } //!< if the returned value is -1, the ged has not been computed

public: /* PUBLIC MEMBER FUNCTIONS */

  /**
   * @brief The method to compute `k` mappings is delegated to derivated classes
   *
   *  This method can include a way to compute a biartite or quadratic cost matrix, to choose
   *  the returned assignments and so on.
   */
  virtual std::list<unsigned int*> getMappings( Graph<NodeAttribute, EdgeAttribute>* g1,
						Graph<NodeAttribute, EdgeAttribute>* g2,
						int k ) = 0;


  /**
   * @brief call to getKOptimalMappings(C, n, m, this->_nep)
   */
  virtual std::list<int*> getKOptimalMappings(double* C,  int n,  int m);

 /**
  * @brief compute an optimal mapping between `g1` and `g2`
  *        from k different optimal mappings by minimizing the ged optained, according
  *        to the ged obtained from `graphdistance`
  * @note The GED is computed and set in `ged`
  */
 virtual double computeOptimalMapping( GraphEditDistance<NodeAttribute,EdgeAttribute> * graphdistance,
                                       Graph<NodeAttribute,EdgeAttribute> * g1,
                                       Graph<NodeAttribute,EdgeAttribute> * g2,
                                       double* C,
                                       unsigned int * G1_to_G2,
				       unsigned int * G2_to_G1 );


 virtual MultiGed<NodeAttribute, EdgeAttribute> * clone() const = 0;

};


// //----


// template<class NodeAttribute, class EdgeAttribute>
// double* MultiGed<NodeAttribute, EdgeAttribute>::
// costMatrixLSAP(double* C, int n, int m)
// {
//   double* Clsap = new  double[(n+m)*(n+m)];
//   //memset(_Clsap, -1.0, sizeof(double)*(n+m)*(n+m)); // inf costs to all non-possible mappings
//   for (int j=0; j<m+n; j++)
//     for (int i=0; i<m+n; i++)
//       if (i>=n && j>=m) Clsap[sub2ind(i,j,n+m)] = 0;
//       else Clsap[sub2ind(i,j,n+m)] = -1.0;

//   //XXX changer column first
//   for(int i=0;i<n;i++)
//     for(int j=0;j<m;j++)
//       Clsap[sub2ind(i,j,n+m)] = C[sub2ind(i,j,n+1)];

//    for(int i=0;i<n;i++)
//      Clsap[sub2ind(i,m+i,n+m)] = C[sub2ind(i,m,n+1)];
//    for(int j=0;j<m;j++)
//      Clsap[sub2ind(n+j,j,n+m)] = C[sub2ind(n,j,n+1)];

//   return Clsap;
// }



template<class NodeAttribute, class EdgeAttribute>
std::list<int*> MultiGed<NodeAttribute, EdgeAttribute>::
getKOptimalMappings( double* C, int n, int m )
{
  return getKOptimalMappings(C, n, m, _nep);
}



template<class NodeAttribute, class EdgeAttribute>
std::list<unsigned int*> MultiGed<NodeAttribute, EdgeAttribute>::
getKOptimalMappings( double* C, int n, int m,  const int& k)
{
  //On considere C fournie
  // compute _Clsap
  // double* Clsap = costMatrixLSAP(C, n, m); //XXX: ? 

  // the returned mappings
  std::list<unsigned int*> solutions;

  __attribute__((unused)) double min_cost = lsape::lsapeSolutions<double>(C,n,m,k,
									  solutions, lsape::ECBP);
  return solutions;
}



template<class NodeAttribute, class EdgeAttribute>
double MultiGed<NodeAttribute, EdgeAttribute>::
computeOptimalMapping ( GraphEditDistance<NodeAttribute,EdgeAttribute> * graphdistance,
                        Graph<NodeAttribute,EdgeAttribute> * g1,
                        Graph<NodeAttribute,EdgeAttribute> * g2,
                        double* C,
                        unsigned int * G1_to_G2, unsigned int * G2_to_G1 )
{

  //XXX:Redondant avec BipartiteGraphEditDistanceLMulti.h
  int n=g1->Size();
  int m=g2->Size();

  std::list<unsigned int*> mappings = getKOptimalMappings(C, n+1, m+1, _nep);

  typename std::list<unsigned int*>::const_iterator it;

  // Get the min of ged;
  _ged = -1;
  double nged;

  unsigned int* local_G1_to_G2 = new unsigned int[n+1];
  unsigned int* local_G2_to_G1 = new unsigned int[m+1];

  for (it=mappings.begin(); it!=mappings.end(); it++){
    unsigned int* lsapMapping = *it;

    nged = graphdistance->GedFromMapping(g1, g2, local_G1_to_G2,n, local_G2_to_G1,m); //XXX: n+1 ?

    // if nged is better : save the mapping and the ged
    if (_ged > nged || _ged == -1){
      _ged = nged;
      for (int i=0; i<n; i++) G1_to_G2[i] = local_G1_to_G2[i]; //XXX:Idem : n+1 ?? quid du dernier élément
      for (int j=0; j<m; j++) G2_to_G1[j] = local_G2_to_G1[j];
    }
  }

  delete [] local_G1_to_G2;
  delete [] local_G2_to_G1;

  typename std::list<unsigned int*>::iterator it_del;
  for (it_del=mappings.begin(); it_del!=mappings.end(); it_del++)
    delete[] *it_del;

  return _ged;
}


#endif

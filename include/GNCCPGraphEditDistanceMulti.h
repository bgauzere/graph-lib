/**
 * @file GNCCPGraphEditDistanceMulti.h
 * @author Évariste <<evariste.daller@unicaen.fr>>
 * @version     0.0.1 -  Mon Apr 10 2017
 *
 */


#ifndef __GNCCPGRAPHEDITDISTANCEMULTI_H__
#define __GNCCPGRAPHEDITDISTANCEMULTI_H__

#include "GNCCPGraphEditDistance.h"
#include "IPFPZetaGraphEditDistanceMulti.h"
#include "MultiGed.h"

template<class NodeAttribute, class EdgeAttribute>
class GNCCPGraphEditDistanceMulti:
  public GNCCPGraphEditDistance<NodeAttribute, EdgeAttribute>,
  public MultiGed<NodeAttribute, EdgeAttribute>
{
protected:

  IPFPZetaGraphEditDistanceMulti<NodeAttribute,EdgeAttribute> * sub_algo;

public: /* CONSTRUCTORS AND ACCESSORS */

  GNCCPGraphEditDistanceMulti( EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction,
                               int n_edit_paths ):
      GNCCPGraphEditDistance<NodeAttribute,EdgeAttribute>(costFunction),
      MultiGed<NodeAttribute,EdgeAttribute>(n_edit_paths)
  {}

  GNCCPGraphEditDistanceMulti( EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction,
                               GraphEditDistance<NodeAttribute,EdgeAttribute> * ed_init,
                               int n_edit_paths ):
    GNCCPGraphEditDistance<NodeAttribute,EdgeAttribute>(costFunction, ed_init),
    MultiGed<NodeAttribute,EdgeAttribute>(n_edit_paths)
  {}


public: /* PUBLIC MEMBER FUNCTIONS */

    virtual void getOptimalMapping( Graph<NodeAttribute,EdgeAttribute> * g1,
                                    Graph<NodeAttribute,EdgeAttribute> * g2,
                                    int * G1_to_G2, int * G2_to_G1 );


};

//----



template<class NodeAttribute, class EdgeAttribute>
void GNCCPGraphEditDistanceMulti<NodeAttribute, EdgeAttribute>::
getOptimalMapping( Graph<NodeAttribute,EdgeAttribute> * g1,
                   Graph<NodeAttribute,EdgeAttribute> * g2,
                   int * G1_to_G2, int * G2_to_G1)
{
  int n = g1->Size();
  int m = g2->Size();

  if(this->_ed_init){
    this->_ed_init->getOptimalMapping(g1,g2,G1_to_G2, G2_to_G1);
  }else{
    for(int i =0;i<g1->Size();i++)
    G1_to_G2[i] = (i>g2->Size())?g2->Size():i;
  for(int j=0;j<g2->Size();j++)
    G2_to_G1[j] = (j>g1->Size())?g1->Size():j;
  }


  this->_zeta = 1;
  this->sub_algo = new IPFPZetaGraphEditDistanceMulti<NodeAttribute, EdgeAttribute>(this->cf, this->_zeta, this->_nep);

  this->sub_algo->setCurrentMatrix(G1_to_G2,G2_to_G1,n,m);
  double * Xk = this->sub_algo->getCurrentMatrix();
  Map<MatrixXd> m_Xk(Xk,n+1,m+1);
#if DEBUG
  IOFormat OctaveFmt(StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
  std::cout << m_Xk.format(OctaveFmt) << std::endl;
#endif
  bool flag = true;
  while((this->_zeta > -1) && flag){
    this->sub_algo->setZeta(this->_zeta);
    this->sub_algo->IPFPalgorithm(g1,g2);
#if DEBUG
    std::cout << "zeta : " << this->_zeta << std::endl;
    std::cout << m_Xk.format(OctaveFmt) << std::endl;
#endif
    this->_zeta -= this->_d;
    flag = (m_Xk.array().round() - m_Xk.array()).abs().sum();
  }

#if DEBUG
  std::cout << "zeta final : " << this->_zeta << std::endl;
  std::cout << m_Xk.format(OctaveFmt) << std::endl;
#endif
  m_Xk= m_Xk * -1;
  this->computeOptimalMapping(this, g1, g2, Xk, G1_to_G2, G2_to_G1);

#if DEBUG
  std::cout << " G1_to_G2 :" << std::endl;
  for (int i = 0; i < n ;i ++)
    std::cout << i << " -> " << G1_to_G2[i] << std::endl;
  std::cout << " G2_to_G1 :" << std::endl;
  for (int j = 0; j < m ;j ++)
    std::cout << j << " -> " << G2_to_G1[j] << std::endl;
#endif
}

#endif //__GNCCPGRAPHEDITDISTANCEMULTI_H__

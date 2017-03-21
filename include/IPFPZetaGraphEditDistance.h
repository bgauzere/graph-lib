/**
 * @file IPFPGraphEditDistance.h
 * @author Benoit <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Wed Mar  8 2017
 * 
 * @todo allow a random init, or alternative initializations
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 */

#ifndef __IPFPZETAGRAPHEDITDISTANCE_H__
#define __IPFPZETAGRAPHEDITDISTANCE_H__
#include <Eigen/Dense>
using namespace Eigen;
#include "hungarian-lsape.hh"
#include "IPFPGraphEditDistance.h"
#include "RandomWalksGraphEditDistance.h"
#include "utils.h"

template<class NodeAttribute, class EdgeAttribute>
class IPFPZetaGraphEditDistance:
  public IPFPGraphEditDistance<NodeAttribute, EdgeAttribute>{
protected:
  
private:
  double _zeta;

public:
  IPFPZetaGraphEditDistance(EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction,
			    GraphEditDistance<NodeAttribute,EdgeAttribute> * ed_init,
			    double zeta):
    IPFPGraphEditDistance<NodeAttribute,EdgeAttribute>(costFunction,ed_init),_zeta(zeta){};
  IPFPZetaGraphEditDistance(EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction,
			    double zeta):
    IPFPGraphEditDistance<NodeAttribute,EdgeAttribute>(costFunction),_zeta(zeta){};
  

  void LinearSubProblem();
  double getCost(int * G1_to_G2,int * G2_to_G1, int n, int m);
  double getAlpha();
  double getBeta();
  
};


template<class NodeAttribute, class EdgeAttribute>
void IPFPZetaGraphEditDistance<NodeAttribute, EdgeAttribute>::
LinearSubProblem(){
  Map<MatrixXd> m_linearSubProblem(this->linearSubProblem,this->n+1,this->m+1);

  Map<MatrixXd> m_XkD(this->XkD,this->n+1,this->m+1);
  Map<MatrixXd> m_C(this->C,this->n+1,this->m+1);
  
  m_linearSubProblem = 2*m_XkD + m_C;
  
}



template<class NodeAttribute, class EdgeAttribute>
double IPFPZetaGraphEditDistance<NodeAttribute, EdgeAttribute>::
getCost(int * G1_to_G2,int * G2_to_G1, int n, int m){
  double S_k = IPFPZetaGraphEditDistance::getCost(G1_to_G2,G2_to_G1, n, m);
  return S_k*(1-fabs(this->zeta)) + this->zeta*linearCost(this->bkp1,this->bkp1,n+1,m+1);
}

template<class NodeAttribute, class EdgeAttribute>
double IPFPZetaGraphEditDistance<NodeAttribute, EdgeAttribute>::
getAlpha(){
  return this->R.back() - 2 * this->S[this->k] + (1-fabs(this->zeta))*this->oldLterm;
}

template<class NodeAttribute, class EdgeAttribute>
double IPFPZetaGraphEditDistance<NodeAttribute, EdgeAttribute>::
getBeta(){
  return this->S.back() + this->S[this->k] - this->R.back() - (1-fabs(this->zeta))*this->oldLterm;
}





#endif // __IPFPZETAGRAPHEDITDISTANCE_H__

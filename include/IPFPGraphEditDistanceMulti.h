/**
 * @file IPFPGraphEditDistance.h
 * @author Evariste <<evariste.daller@unicaen.fr>>
 * @version     0.0.1 - Mon Apr  3 2017
 *
 *
 * Description of the program objectives.
 * All necessary references.
 */

#ifndef __IPFPGRAPHEDITDISTANCEMULTI_H__
#define __IPFPGRAPHEDITDISTANCEMULTI_H__

#include <list>
#include "IPFPGraphEditDistance.h"
#include "MultiGed.h"


using namespace Eigen;

template<class NodeAttribute, class EdgeAttribute>
class IPFPGraphEditDistanceMulti:
  public virtual IPFPGraphEditDistance<NodeAttribute, EdgeAttribute>,
  public MultiGed<NodeAttribute, EdgeAttribute>
{

public:

  virtual void IPFPalgorithm( Graph<NodeAttribute,EdgeAttribute> * g1,
                              Graph<NodeAttribute,EdgeAttribute> * g2 );


  virtual void getOptimalMapping( Graph<NodeAttribute,EdgeAttribute> * g1,
				                          Graph<NodeAttribute,EdgeAttribute> * g2,
				                          int * G1_to_G2, int * G2_to_G1 );



  IPFPGraphEditDistanceMulti( EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction,
                              GraphEditDistance<NodeAttribute,EdgeAttribute> * ed_init,
                              int n_edit_paths
                            ):
    IPFPGraphEditDistance<NodeAttribute,EdgeAttribute>(costFunction, ed_init),
    MultiGed<NodeAttribute,EdgeAttribute>(n_edit_paths)
  {}


};

//---


template<class NodeAttribute, class EdgeAttribute>
void IPFPGraphEditDistanceMulti<NodeAttribute, EdgeAttribute>::
getOptimalMapping( Graph<NodeAttribute,EdgeAttribute> * g1,
		               Graph<NodeAttribute,EdgeAttribute> * g2,
		               int * G1_to_G2, int * G2_to_G1)
{
  //Compute Mapping init
  if (this->_ed_init)
    this->_ed_init->getOptimalMapping(g1,g2,G1_to_G2, G2_to_G1);
  this->_n = g1->Size();
  this->_m = g2->Size();

  this->Xk = new double[(this->_n+1)*(  this->_m+1)];
  Map<MatrixXd> m_Xk(this->Xk,  this->_n+1,  this->_m+1);
  this->Xk = this->mappingsToMatrix(G1_to_G2, G2_to_G1, this->_n, this->_m, this->Xk);

  this->IPFPalgorithm(g1, g2);


  m_Xk = m_Xk *-1;
  this->computeOptimalMapping(this, g1, g2, this->Xk, G1_to_G2, G2_to_G1);
  delete [] this->Xk; this->Xk = NULL;
}



template<class NodeAttribute, class EdgeAttribute>
void IPFPGraphEditDistanceMulti<NodeAttribute, EdgeAttribute>::
IPFPalgorithm(Graph<NodeAttribute,EdgeAttribute> * g1,
              Graph<NodeAttribute,EdgeAttribute> * g2){

  this->_directed = (g1->isDirected() && g2->isDirected());
  //We assume that Xk is filled with a matrix, binary or not

  this->_n = g1->Size();
  this->_m = g2->Size();

  this->S.clear();
  this->R.clear();

  this->NodeCostMatrix(g1,g2);//REdondant for GNCCP
  Map<MatrixXd> m_C(this->C,this->_n+1,this->_m+1); //REdondant for GNCCP

  this->bkp1 = new double [(  this->_n+1) * (  this->_m+1)];
  Map<MatrixXd> m_bkp1 (this->bkp1,  this->_n+1,  this->_m+1);
  //this->bkp1 = mappingsToMatrix(G1_to_G2,G2_to_G1,  this->_n,  this->_m,this->bkp1);

  this->XkD = this->QuadraticTerm(g1,g2,this->Xk,NULL); //REdondant for GNCCP
  Map<MatrixXd> m_XkD(this->XkD,  this->_n+1,  this->_m+1);

  Map<MatrixXd> m_Xk(this->Xk,  this->_n+1,  this->_m+1);

  this->Lterm = this->linearCost(this->C,this->Xk,  this->_n+1,  this->_m+1);
  this->S.push_back(this->getCost(this->Xk,  this->_n,  this->_m));
#if DEBUG
  std::cout << "S(0) = " << this->S.back() << std::endl;
#endif
  this->k=0;
  this->linearSubProblem = new double [(  this->_n+1) * (  this->_m+1)];
  Map<MatrixXd> m_linearSubProblem(this->linearSubProblem,  this->_n+1,  this->_m+1);


  this->Xkp1tD = new double [(  this->_n+1) * (  this->_m+1)];
  Map<MatrixXd> m_Xkp1tD (this->Xkp1tD,  this->_n+1,  this->_m+1);


  //double *u = new double[this->_n+1];
  //double *v = new double[this->_m+1];
  int * G1_to_G2 = new int[this->_n];
  int * G2_to_G1 = new int[this->_m];
  bool flag_continue = true;

  while((this->k < this->maxIter) && flag_continue){ //TODO : fixer un epsilon, param ?
    this->XkD = this->QuadraticTerm(g1,g2,this->Xk,this->XkD);
    this->LinearSubProblem();//    should call it gradient direction

    //hungarianLSAPE(linearSubProblem,  this->_n+1,  this->_m+1, G1_to_G2,G2_to_G1, u,v,false);
    this->computeOptimalMapping(this, g1, g2, this->linearSubProblem, G1_to_G2, G2_to_G1);

    //bkp1 is the matrix version of mapping G1_to_G2 and G2_to_G1, so a binary matrix
    this->bkp1 = this->mappingsToMatrix(G1_to_G2,G2_to_G1,  this->_n,  this->_m,this->bkp1);
    this->R.push_back(this->linearCost(this->linearSubProblem,G1_to_G2, G2_to_G1,  this->_n,  this->_m));

#if DEBUG
    IOFormat OctaveFmt(StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
    std::cout << "XkD" << std::endl;
    std::cout << m_XkD.format(OctaveFmt) << std::endl;
    std::cout << "linearSubProblem" << std::endl;
    std::cout << m_linearSubProblem.format(OctaveFmt) << std::endl;
    std::cout << "bkp1" << std::endl;
    std::cout << m_bkp1.format(OctaveFmt) << std::endl;
    std::cout << "R : " << this->R.back() << std::endl;
#endif

    this->oldLterm = this->Lterm;
    this->Lterm = this->linearCost(this->C,G1_to_G2, G2_to_G1,  this->_n,  this->_m);
    this->XkD = this->QuadraticTerm(g1,g2,G1_to_G2, G2_to_G1,this->XkD);
    this->S.push_back(this->getCost(G1_to_G2, G2_to_G1,  this->_n,  this->_m));

#if DEBUG
    std::cout << "S : " << this->S.back() << std::endl;

#endif

    double alpha = this->getAlpha();
    double beta= this->getBeta();
    double t0 =0.0;
    if(beta > 0.000001)
      t0 = -alpha / (2.*beta);
    //Built a new Xk matrix (possibly not permutation)
#if DEBUG

      std::cout << "t0 : " << t0 << std::endl;
      std::cout << "Alpha : " << alpha << std::endl;
      std::cout << "Beta : " << beta << std::endl;
#endif

    if ((beta < 0.00001) || (t0 >= 1)){
      //Check if Xk evolves
      if((m_Xk - m_bkp1).norm() < 0.0001){
	flag_continue = false;
      }else{
	memcpy(this->Xk, this->bkp1,sizeof(double)*(  this->_n+1)*(  this->_m+1));
      }
      //Lterm = Lterm_new;
    }else{
      //Line search
      MatrixXd maj_matrix(  this->_n+1,  this->_m+1);
      maj_matrix = t0*(m_bkp1 - m_Xk);
#if DEBUG
      std::cout << "line search" << std::endl;
      std::cout << "Norm de la maj : " << maj_matrix.norm() << std::endl;
#endif
      if(maj_matrix.norm() < 0.0001){
	flag_continue = false;
      }else {
	m_Xk = m_Xk + t0*(m_bkp1 - m_Xk);
	this->S[this->k+1] = this->S[this->k] - ((pow(alpha,2))/(4*beta));
	this->Lterm = this->linearCost(this->C, this->Xk,   this->_n+1,  this->_m+1);
      }
    }
#if DEBUG
    std::cout << "Xk à l'itération " << this->k << std::endl;
    std::cout << m_Xk.format(OctaveFmt) << std::endl;
    std::cout << "------------------------------------------------------------"  << std::endl << std::endl;
#endif

    this->k++;
  }
  delete [] this->Xkp1tD;this->Xkp1tD=0;
  delete [] this->linearSubProblem;this->linearSubProblem=0;
  delete [] this->XkD;this->XkD = 0;
  delete [] this->C;this->C = 0;
  delete [] this->bkp1;this->bkp1 = 0;
  //delete [] u;
  //delete [] v;
  delete [] G1_to_G2;
  delete [] G2_to_G1;
#if DEBUG
  std::cout << this->S.back() << std::endl;
  std::cout << "Fin d'IPFP : "<< this->k << "iterations " << std::endl;
#endif
  //Xk contains the optimal bistochastic matrix, binary or not.

}


#endif // __IPFPGRAPHEDITDISTANCEMULTI_H__

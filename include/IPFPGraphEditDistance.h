/**
 * @file IPFPGraphEditDistance.h
 * @author Benoit <<benoit.gauzere@insa-rouen.fr>>
 * @version     0.0.1 - Wed Mar  8 2017
 *
 * @todo allow a random init, or alternative initializations
 * @bug the list of known bugs.
 *
 * Description of the program objectives.
 * All necessary references.
 */

#ifndef __IPFPGRAPHEDITDISTANCE_H__
#define __IPFPGRAPHEDITDISTANCE_H__
#include <Eigen/Dense>
#include <limits>
using namespace Eigen;

#include <lsape.h> // Bistochastic generation and sinkhorn balancing
#include "GraphEditDistance.h"
#include "MultiGed.h"
#include "IPFPQAP.h"
#include "gl_utils.h"

template<class NodeAttribute, class EdgeAttribute>
class IPFPGraphEditDistance:
  public IPFPQAP<NodeAttribute, EdgeAttribute>,
  public GraphEditDistance<NodeAttribute, EdgeAttribute>
{

protected:

  GraphEditDistance<NodeAttribute,EdgeAttribute> * _ed_init;
  bool cleanCostFunction;
  bool useContinuousRandomInit;
  bool useContinuousFlatInit;
  bool useSinkhorn;

  virtual
  void NodeCostMatrix(Graph<NodeAttribute,EdgeAttribute> * g1,
		      Graph<NodeAttribute,EdgeAttribute> * g2);

  virtual
  double * QuadraticTerm(Graph<NodeAttribute,EdgeAttribute> * g1,
			 Graph<NodeAttribute,EdgeAttribute> * g2,
			 unsigned int * G1_to_G2, unsigned int * G2_to_G1, double * XkD);

  virtual
  double * QuadraticTerm(Graph<NodeAttribute,EdgeAttribute> * g1,
                         Graph<NodeAttribute,EdgeAttribute> * g2,
                         double * Matrix, double * XkD);
  virtual
  double * QuadraticTerm(Graph<NodeAttribute,EdgeAttribute> * g1,
			 Graph<NodeAttribute,EdgeAttribute> * g2,
			 std::vector<std::pair<std::pair<unsigned int,unsigned int>,double> > mappings, double * XkD);

  //This linearCost is efficient for sparse Xk matrices
  // n is nb rows of matrices, m is nb columns
  virtual double linearCost(double * CostMatrix, unsigned int * G1_to_G2,unsigned int * G2_to_G1, unsigned int n, unsigned int m);
  // n is nb rows of matrices, m is nb columns
  virtual double linearCost(double * CostMatrix, double * Xk, unsigned int n, unsigned int m);

  // Fill this->linearSubProblem with appropriatelinear problem
  virtual void LinearSubProblem();
  virtual double getCost(unsigned int * G1_to_G2,unsigned int * G2_to_G1, unsigned int n, unsigned int m);
  virtual double getCost(double * Matrix , unsigned int n, unsigned int m);
  virtual double getAlpha();
  virtual double getBeta();

  virtual double * mappingsToMatrix(unsigned int * G1_to_G2,unsigned int * G2_to_G1, unsigned int n, unsigned int m, double * Matrix);



public:
  IPFPGraphEditDistance(EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction,
			GraphEditDistance<NodeAttribute,EdgeAttribute> * ed_init):
    IPFPQAP<NodeAttribute, EdgeAttribute>(costFunction),
    GraphEditDistance<NodeAttribute,EdgeAttribute>(costFunction),
    _ed_init(ed_init),
    cleanCostFunction(false),
    useContinuousRandomInit(false),
    useContinuousFlatInit(false),
    useSinkhorn(false)
  {};

  IPFPGraphEditDistance(EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction):
    IPFPQAP<NodeAttribute, EdgeAttribute>(costFunction),
    GraphEditDistance<NodeAttribute,EdgeAttribute>(costFunction),
    _ed_init(NULL),
    cleanCostFunction(false),
    useContinuousRandomInit(false),
    useContinuousFlatInit(false),
    useSinkhorn(false)
  {
    this->C = NULL; this->linearSubProblem=NULL; this->XkD=NULL; this->Xk=NULL; this->Lterm=0; this->oldLterm=0;
    this->Xkp1tD=NULL; this->bkp1=NULL; this->_n=-1; this->_m=-1; this->k=-1; this->_directed=false;
    this->_compute_equiv_mappings = false;
  };


 /**
  * @brief  Compute in `G1_to_G2` and `G2_to_G1` the best mapping that the
  *         method found between `g1` and `g2`
  */
  virtual void getOptimalMapping(Graph<NodeAttribute,EdgeAttribute> * g1,
				 Graph<NodeAttribute,EdgeAttribute> * g2,
				 unsigned int * G1_to_G2, unsigned int * G2_to_G1);


  /**
   * @brief  Refine the given mapping in `G1_to_G2`, `G2_to_G1` between
   *         the graphs `g1` and `g2`. The result replaces the initial
   *         values of these arrays.
   *
   *  If the flag _compute_equiv_mappings is true, the two arrays
   *  `G1_to_G2` and `G2_to_G1` keep their initial values and
   *  the refined mappings can be read through the methods
   *  getEquivalentG2toG1 and getEquivalentG1toG2.
   *
   * @param  fromInit  Ignored
   */
  virtual void getBetterMapping( Graph<NodeAttribute,EdgeAttribute> * g1,
				 Graph<NodeAttribute,EdgeAttribute> * g2,
				 unsigned int * G1_to_G2, unsigned int * G2_to_G1, bool fromInit=true);

  /**
   * @brief  Perform the Frank-Wlofe algorithm
   */
  void IPFPalgorithm(Graph<NodeAttribute,EdgeAttribute> * g1,
		     Graph<NodeAttribute,EdgeAttribute> * g2);



  virtual double mappingCost( Graph<NodeAttribute,EdgeAttribute> * g1,
                      Graph<NodeAttribute,EdgeAttribute> * g2,
                     unsigned int* G1_to_G2, unsigned int* G2_to_G1 );


  /**
   * @brief  Activate the use of a continuous bistochastic random matrix as initialization, no matter what
   *         is provided to the `get****Mapping()` method.
   */
  void continuousRandomInit(bool yes=true){
    this->useContinuousRandomInit = yes;
  }


  /**
   * @brief  Activate or deactivate the continuous uniform flat initialization, no matter what is provided.
   *
   *   The geometrical barycenter of error correcting doubly stochastic matrices
   *   \f$J = \frac{2(n+1)\times(m+1)}{n+m+2}\f$ will be used as initialization in any case
   */
  void continuousFlatInit(bool yes=true){
    this->useContinuousFlatInit = yes;
  }

  /*
   * @brief  Activate the use of a Sinkhorn balancing of the initialization before starting IPFP algorithm
   *
   * @note  If the flags `recenter` and `useSinjhorn` are both `true`,
   *        then the Sinkhorn balancing is done before the centering procedure.
   *
  void sinkhornFirst(bool yes=true){
    this->useSinkhorn = yes;
  }
  //*/

  void recenterInit(bool yes=true){
    this->recenter = yes;
  }

  /**
   * @brief  Set the centering matrix to J of size \f$(n+1)\times(m+1)\f$ and activate the centering
   *
   *  On the next IPFP algorithm call, the initialization \f$X_0\f$ will be translated to \f$\frac{1}{2}\times(X_0+J)\f$.
   *  If `nJ==NULL`, then the geometrical barycenter of error correcting doubly stochastic matrices
   *   \f$J = \frac{2(n+1)\times(m+1)}{n+m+2}\f$ will be used
   */
  virtual void recenterInit(double* nJ, unsigned int n, unsigned int m);

  IPFPGraphEditDistance * clone() const { return new IPFPGraphEditDistance(*this); }

  IPFPGraphEditDistance( const IPFPGraphEditDistance& other ) :
    IPFPQAP<NodeAttribute,EdgeAttribute>(NULL),
    GraphEditDistance<NodeAttribute, EdgeAttribute>(NULL),
    _ed_init(other._ed_init),
    cleanCostFunction(true)
  {
    this->cf = other.cf->clone();
    this->costFunction = this->cf;
    this->recenter = other.recenter;
    this->useContinuousRandomInit = other.useContinuousRandomInit;
    this->useContinuousFlatInit = other.useContinuousFlatInit;
    this->useSinkhorn = other.useSinkhorn;
    this->J=NULL;
    this->_compute_equiv_mappings = other._compute_equiv_mappings;
  }

  virtual ~IPFPGraphEditDistance(){
    if (this->C != NULL) delete [] this->C;
    if (this->linearSubProblem != NULL) delete [] this->linearSubProblem;
    if (this->XkD != NULL) delete [] this->XkD;
    if (this->Xk != NULL) delete [] this->Xk;
    if (this->Xkp1tD != NULL) delete [] this->Xkp1tD;
    if (this->bkp1 != NULL) delete [] this->bkp1;
    if (this->cleanCostFunction) delete this->cf;
  }

};


template<class NodeAttribute, class EdgeAttribute>
void IPFPGraphEditDistance<NodeAttribute, EdgeAttribute>::
recenterInit(double * nJ, unsigned int n, unsigned int m)
{
  if ( this->J != NULL) delete[] this->J;
  this->J = new double[(n+1)*(m+1)];

  this->recenter = true;
  if (nJ == NULL){
    for (unsigned int j=0; j<m+1; j++)
      for (unsigned int i=0; i<n+1; i++)
        this->J[sub2ind(i,j,n+1)] = 2.0 / (n+m+2);
        // this->J[sub2ind(i,j,n+1)] = 1.0 / ((n+1)*(m+1) -1) ;
  }
  else{
    for (unsigned int j=0; j<m+1; j++)
      for (unsigned int i=0; i<n+1; i++)
        this->J[sub2ind(i,j,n+1)] = nJ[sub2ind(i,j,n+1)];
  }
}


template<class NodeAttribute, class EdgeAttribute>
void IPFPGraphEditDistance<NodeAttribute, EdgeAttribute>::NodeCostMatrix(Graph<NodeAttribute,EdgeAttribute> * g1,
									 Graph<NodeAttribute,EdgeAttribute> * g2){
  unsigned int n=g1->Size();
  unsigned int m=g2->Size();

  this->C = new double[(n+1)*(m+1)];
  //memset(this->C,std::numeric_limits<double>::max(),sizeof(double)*(n+1)*(m+1));
  this->C[sub2ind(n,m,(n+1))] = 0;
  for(unsigned int i=0;i<n;i++)
    for(unsigned int j=0;j<m;j++)
      this->C[sub2ind(i,j,(n+1))] = this->cf->NodeSubstitutionCost((*g1)[i],(*g2)[j],g1,g2);

  for(unsigned int i=0;i<n;i++)
    this->C[sub2ind(i,m,(n+1))] = this->cf->NodeDeletionCost((*g1)[i],g1);

  for(unsigned int j=0;j<m;j++)
    this->C[sub2ind(n,j,(n+1))] = this->cf->NodeInsertionCost((*g2)[j],g2);
}


template<class NodeAttribute, class EdgeAttribute>
double * IPFPGraphEditDistance<NodeAttribute,
			       EdgeAttribute>::QuadraticTerm(Graph<NodeAttribute,EdgeAttribute> * g1,
							     Graph<NodeAttribute,EdgeAttribute> * g2,
							     double * Matrix, double * XkD){
  unsigned int n = g1->Size();
  unsigned int m = g2->Size();

  std::vector<std::pair<std::pair<unsigned int,unsigned int>, double>> mappings;
  for(unsigned int i=0;i<n+1;i++)
    for(unsigned int j=0;j<m+1;j++){
      double value = Matrix[sub2ind(i,j,n+1)];
      if(value > 0.){
	std::pair<unsigned int,unsigned int> tmp = std::pair<unsigned int,unsigned int>(i,j);
	mappings.push_back(std::pair<std::pair<unsigned int,unsigned int>,double>(tmp,value));
      }
    }
  if (!(this->XkD))
    this->XkD=new double[(n+1)*(m+1)];

  return this->QuadraticTerm(g1,g2,mappings, this->XkD);

}



template<class NodeAttribute, class EdgeAttribute>
double * IPFPGraphEditDistance<NodeAttribute,
			       EdgeAttribute>::QuadraticTerm(Graph<NodeAttribute,EdgeAttribute> * g1,
							     Graph<NodeAttribute,EdgeAttribute> * g2,
							     unsigned int * G1_to_G2, unsigned int * G2_to_G1,double * XkD){

  unsigned int n = g1->Size();
  unsigned int m = g2->Size();

  //Reconstruction d'un mapping
  std::vector<std::pair<std::pair<unsigned int,unsigned int>, double>> mappings;
  for (unsigned int i =0;i<n;i++){
    std::pair<int,int> tmp = std::pair<unsigned int,unsigned int>(i,G1_to_G2[i]);
    mappings.push_back(std::pair<std::pair<unsigned int,unsigned int>,double>(tmp,1.));
  }
  for (unsigned int j =0;j<m;j++)
    if (G2_to_G1[j] >= n){
      std::pair<int,int> tmp = std::pair<unsigned int,unsigned int>(G2_to_G1[j],j);
      mappings.push_back(std::pair<std::pair<unsigned int,unsigned int>,double>(tmp,1.));
    }
  if (!(this->XkD))
    XkD=new double[(n+1)*(m+1)];

  return this->QuadraticTerm(g1,g2,mappings,this->XkD);

}



template<class NodeAttribute, class EdgeAttribute>
double * IPFPGraphEditDistance<NodeAttribute,
			       EdgeAttribute>::QuadraticTerm(Graph<NodeAttribute,EdgeAttribute> * g1,
							     Graph<NodeAttribute,EdgeAttribute> * g2,
							     std::vector<std::pair<std::pair<unsigned int,unsigned int>,double> > mappings,
							     double * quadraticTerm){
  unsigned int n = g1->Size();
  unsigned int m = g2->Size();

  if (! quadraticTerm)
    quadraticTerm=new double[(n+1)*(m+1)];

  memset(quadraticTerm,0,sizeof(double)*(n+1)*(m+1));

  for(unsigned int j = 0; j < n+1; j++){ // Attention : dans le papier sspr, condition sur x_jl /= 0. En effet, inutile pour le cas ou on multiplie a droite par le mapping. Mais nécessaire quand on utilise XtD dans le sous probleme
    for(unsigned int l = 0; l < m+1;l++){

      std::vector<std::pair<std::pair<unsigned int,unsigned int>,double> >::iterator it = mappings.begin();
      for(;it != mappings.end();it++){
	unsigned int i = it->first.first;
	unsigned int k = it->first.second;
	bool eps_i,eps_j,eps_k,eps_l;
	eps_i = (i >= n);eps_j = (j >= n);eps_k = (k >= m);eps_l = (l >= m);

	if( ((i != j) || eps_i) && ((k != l) || eps_k)){
	  GEdge<EdgeAttribute> * e1 = NULL;
	  bool delta_e1 = false;
	  if ((!eps_i) && (!eps_j)){
	    e1 = g1->getEdge(i,j);
	    delta_e1 = (e1 !=NULL);
	  }

	  GEdge<EdgeAttribute> * e2 = NULL;
	  bool delta_e2 = false;
	  if((! eps_k) && (! eps_l)){
	    e2=g2->getEdge(k,l);
	    delta_e2 = (e2 != NULL);// false if l>m
	  }
	  //TODO : Optimize if sequence
	  //If (i,j) and (k,l) are both same nodes,
	  //no edges between them, so delta_e1 and delta_e2 are both 0, and so the cost

	  double cost = 0.0;
	  if (delta_e1 && delta_e2) // sub
	    cost = this->cf->EdgeSubstitutionCost(e1,e2,g1,g2);
	  else if ((delta_e1) && (!delta_e2)) //deletion
	    cost = this->cf->EdgeDeletionCost(e1,g1);
	  else if ((! delta_e1) && delta_e2)
	    cost = this->cf->EdgeInsertionCost(e2,g2);

	  quadraticTerm[sub2ind(j,l,n+1)] +=  cost*it->second;

	}
// #if DEBUG
// 	std::cout << "i : " << i<< std::endl;
// 	std::cout << "j : " << j<< std::endl;
// 	std::cout << "k : " << k<< std::endl;
// 	std::cout << "l : " << l<< std::endl;

// 	std::cout << "eps i : " << eps_i<< std::endl;
// 	std::cout << "eps j : " << eps_j<< std::endl;
// 	std::cout << "eps k : " << eps_k<< std::endl;
// 	std::cout << "eps l : " << eps_l<< std::endl;

// 	std::cout << "delta_e1 : " << delta_e1<< std::endl;
// 	std::cout << "delta_e2 : " << delta_e2<< std::endl;

// 	std::cout << "cost : " << cost << std::endl;

// #endif
      }
      if(! this->_directed)
      	quadraticTerm[sub2ind(j,l,n+1)] *= 0.5;

    }
  }
  return quadraticTerm;


}


template<class NodeAttribute, class EdgeAttribute>
void IPFPGraphEditDistance<NodeAttribute, EdgeAttribute>::
getOptimalMapping( Graph<NodeAttribute,EdgeAttribute> * g1,
                   Graph<NodeAttribute,EdgeAttribute> * g2,
                   unsigned int * G1_to_G2, unsigned int * G2_to_G1 )
{
  //Compute Mapping init
  if (this->_ed_init)
    this->_ed_init->getOptimalMapping(g1,g2,G1_to_G2, G2_to_G1);

  getBetterMapping(g1, g2, G1_to_G2, G2_to_G1);
}



template<class NodeAttribute, class EdgeAttribute>
void IPFPGraphEditDistance<NodeAttribute, EdgeAttribute>::
getBetterMapping( Graph<NodeAttribute,EdgeAttribute> * g1,
                  Graph<NodeAttribute,EdgeAttribute> * g2,
                  unsigned int * G1_to_G2, unsigned int * G2_to_G1, bool fromInit)
{
  this->_n = g1->Size();
  this->_m = g2->Size();

  this->Xk = new double[(this->_n+1)*(  this->_m+1)];
  Map<MatrixXd> m_Xk(this->Xk,  this->_n+1,  this->_m+1);

  if (!useContinuousRandomInit && !useContinuousFlatInit)
    this->Xk = this->mappingsToMatrix(G1_to_G2,G2_to_G1,this->_n,this->_m,this->Xk);

  this->IPFPalgorithm(g1,g2);


  m_Xk= MatrixXd::Ones(this->_n+1, this->_m+1)-m_Xk;

  if (this->_compute_equiv_mappings){
    // Get less precision to induce more equivalent mappings
    m_Xk *= 100;
    for (unsigned int i=0; i<(this->_n+1)*(this->_m+1); i++) this->Xk[i] = (int)(this->Xk[i]);

    //TODO virer le 100
    std::list<unsigned int*> lsapEquivMaps = MultiGed<NodeAttribute, EdgeAttribute>::getKOptimalMappings(this->Xk, this->_n, this->_m, 100);
    MultiGed<NodeAttribute, EdgeAttribute>::mappings_lsap2lsape(
      lsapEquivMaps,
      this->_equiv_G1toG2, this->_equiv_G2toG1,
      this->_n, this->_m
    );

    for (std::list<unsigned int*>::iterator it=lsapEquivMaps.begin(); it != lsapEquivMaps.end();  it++)
      delete [] (*it);
  }
  else{
    double *u = new double[this->_n+1];
    double *v = new double[this->_m+1];
    lsape::hungarianLSAPE(this->Xk,  this->_n+1,  this->_m+1, G1_to_G2,G2_to_G1, u,v,false);

    // Save also the mapping in the list of equivalent mappings
    this->_equiv_G1toG2.clear();
    this->_equiv_G2toG1.clear();
    this->_equiv_G1toG2.push_back(G1_to_G2);
    this->_equiv_G2toG1.push_back(G2_to_G1);
    delete [] u;
    delete [] v;
  }
  delete [] this->Xk; this->Xk=0;
}

template<class NodeAttribute, class EdgeAttribute>
void IPFPGraphEditDistance<NodeAttribute, EdgeAttribute>::
IPFPalgorithm(Graph<NodeAttribute,EdgeAttribute> * g1,
	      Graph<NodeAttribute,EdgeAttribute> * g2)
{

  // If use a random bistochastic continuous matrix :
  // MODIF
  // if (useContinuousRandomInit){
  //   double* _I = randBiStochExt<double,int>(this->_n, this->_m);
  //   reduceExt(_I, this->_n, this->_m, this->Xk);
  //   delete [] _I;
  // }

  // Use J as init
  if (useContinuousFlatInit){
    this->recenterInit(NULL, this->_n, this->_m);
    this->recenter = false;
    for (unsigned int j=0; j<this->_m+1;  j++)
      for (unsigned int i=0; i<this->_n+1;  i++)
        this->Xk[sub2ind(i,j,this->_n+1)] = this->J[sub2ind(i,j,this->_n+1)];
  }

  // Recenter the init mapping ?
  if (this->recenter){
    if (this->J == NULL) this->recenterInit(NULL, this->_n, this->_m);
    for (unsigned int j=0; j<this->_m+1;  j++)
      for (unsigned int i=0; i<this->_n+1;  i++)
        this->Xk[sub2ind(i,j,this->_n+1)] = (this->Xk[sub2ind(i,j,this->_n+1)] + this->J[sub2ind(i,j,this->_n+1)]) / 2;
  }

  this->_directed = (g1->isDirected() && g2->isDirected());
  //We assume that Xk is filled with a matrix, binary or not

  this->_n = g1->Size();
  this->_m = g2->Size();

  this->S.clear();
  this->R.clear();

  NodeCostMatrix(g1,g2);//REdondant for GNCCP
  Map<MatrixXd> m_C(this->C,this->_n+1,this->_m+1); //REdondant for GNCCP

  this->bkp1 = new double [(  this->_n+1) * (  this->_m+1)];
  Map<MatrixXd> m_bkp1 (this->bkp1,  this->_n+1,  this->_m+1);
  //this->bkp1 = mappingsToMatrix(G1_to_G2,G2_to_G1,  this->_n,  this->_m,this->bkp1);

  this->XkD = this->QuadraticTerm(g1,g2,this->Xk,NULL); //REdondant for GNCCP
  Map<MatrixXd> m_XkD(this->XkD,  this->_n+1,  this->_m+1);

  Map<MatrixXd> m_Xk(this->Xk,  this->_n+1,  this->_m+1);

  this->Lterm = linearCost(this->C,this->Xk,  this->_n+1,  this->_m+1);
  this->S.push_back(this->getCost(this->Xk,  this->_n,  this->_m));
#if DEBUG
  std::cout << "S(0) = " << this->S.back() << std::endl;
#endif
  this->k=0;
  this->linearSubProblem = new double [(  this->_n+1) * (  this->_m+1)];
  Map<MatrixXd> m_linearSubProblem(this->linearSubProblem,  this->_n+1,  this->_m+1);


  this->Xkp1tD = new double [(  this->_n+1) * (  this->_m+1)];
  Map<MatrixXd> m_Xkp1tD (this->Xkp1tD,  this->_n+1,  this->_m+1);


  double *u = new double[this->_n+1];
  double *v = new double[this->_m+1];  unsigned int * G1_to_G2 = new unsigned int[this->_n];
  unsigned int * G2_to_G1 = new unsigned int[this->_m];
  bool flag_continue = true;

  //BipartiteGraphEditDistanceMulti<int,int> ed_multi(this->cf, 30); // To know how many solutions to lsap per iteration
  while((this->k < this->maxIter) && flag_continue){ //TODO : fixer un epsilon, param ?
    this->XkD = QuadraticTerm(g1,g2,this->Xk,this->XkD);
    this->LinearSubProblem();//    should call it gradient direction

    lsape::hungarianLSAPE(this->linearSubProblem,  this->_n+1,  this->_m+1, G1_to_G2,G2_to_G1, u,v,false);
    //bkp1 is the matrix version of mapping G1_to_G2 and G2_to_G1, so a binary matrix
    this->bkp1 = mappingsToMatrix(G1_to_G2,G2_to_G1,  this->_n,  this->_m,this->bkp1);
    this->R.push_back(linearCost(this->linearSubProblem,G1_to_G2, G2_to_G1,  this->_n,  this->_m));
    //std::list<int*> mappings = ed_multi.getKOptimalMappings(g1, g2, linearSubProblem, 30);
    //std::cout << mappings.size() << ", ";

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
    this->Lterm = linearCost(this->C,G1_to_G2, G2_to_G1,  this->_n,  this->_m);
    this->XkD = QuadraticTerm(g1,g2,G1_to_G2, G2_to_G1, this->XkD);
    this->S.push_back(this->getCost(G1_to_G2, G2_to_G1,  this->_n,  this->_m));

#if DEBUG
    std::cout << "S : " << this->S.back() << std::endl;

#endif

    double alpha = getAlpha();
    double beta= getBeta();
    double t0 =0.0;
    if(beta > 0.000001)
      t0 = -alpha / (2.*beta);
    //Built a new Xk matrix (possibly not permutation)
#if DEBUG

      std::cout << "t0 : " << t0 << std::endl;
      std::cout << "Alpha : " << alpha << std::endl;
      std::cout << "Beta : " << beta << std::endl;
#endif

    /*
    if (R.back() < 0.0001 || -alpha / R.back() < this->epsilon){
      flag_continue = false;
    }
    //*/

    //*
    if (this->R.back() < 0.0001)
      flag_continue = (fabs(alpha) > this->epsilon);
    else
      flag_continue = (fabs(alpha / this->R.back()) > this->epsilon);
    //*/

    if ((beta < 0.00001) || (t0 >= 1))
      //if(flag_continue)
        memcpy(this->Xk, this->bkp1,sizeof(double)*(  this->_n+1)*(  this->_m+1));

      //Lterm = Lterm_new;
    else{
      //Line search
      MatrixXd maj_matrix(  this->_n+1,  this->_m+1);
      maj_matrix = t0*(m_bkp1 - m_Xk);
#if DEBUG
      std::cout << "line search" << std::endl;
      std::cout << "Norm de la maj : " << maj_matrix.norm() << std::endl;
#endif
      //if (flag_continue){
        m_Xk = m_Xk + t0*(m_bkp1 - m_Xk);
        this->S[this->k+1] = this->S[this->k] - ((pow(alpha,2))/(4*beta));
        this->Lterm = linearCost(this->C, this->Xk,   this->_n+1,  this->_m+1);
	//}
    }
#if DEBUG
    std::cout << "Xk à l'itération " << k << std::endl;
    std::cout << m_Xk.format(OctaveFmt) << std::endl;
    std::cout << "------------------------------------------------------------"  << std::endl << std::endl;
#endif

    this->k++;
  }

  //std::cout << this->k << ", " ;

  delete [] this->Xkp1tD;this->Xkp1tD=0;
  delete [] this->linearSubProblem;this->linearSubProblem=0;
  delete [] this->XkD;this->XkD = 0;
  delete [] this->C;this->C = 0;
  delete [] this->bkp1;this->bkp1 = 0;
  delete [] u;
  delete [] v;
  delete [] G1_to_G2;
  delete [] G2_to_G1;
#if DEBUG
  std::cout << this->S.back() << std::endl;
  std::cout << "Fin d'IPFP : "<< k << "iterations " << std::endl;
#endif
  //Xk contains the optimal bistochastic matrix, binary or not.

}



template<class NodeAttribute, class EdgeAttribute>
double IPFPGraphEditDistance<NodeAttribute, EdgeAttribute>::
linearCost(double * CostMatrix, unsigned int * G1_to_G2,unsigned int * G2_to_G1, unsigned int n, unsigned int m){
  double sum = 0.0;
  for(unsigned int i=0;i<n;i++)
    sum += CostMatrix[sub2ind(i,G1_to_G2[i],n+1)];
  for(unsigned int j=0;j<m;j++)
    if(G2_to_G1[j] >= n)
      sum+= CostMatrix[sub2ind(G2_to_G1[j],j,n+1)];
  return sum;
}

template<class NodeAttribute, class EdgeAttribute>
double IPFPGraphEditDistance<NodeAttribute, EdgeAttribute>::
linearCost(double * CostMatrix, double * X, unsigned int n, unsigned int m){
  //Todo : optimiser avec dot product ?
  double sum = 0.0;
  for(unsigned int i=0;i<n;i++)
    for(unsigned int j=0;j<m;j++)
      sum += CostMatrix[sub2ind(i,j,n)] * X[sub2ind(i,j,n)]; //XXX is that not n+1 ?
  return sum;
}

template<class NodeAttribute, class EdgeAttribute>
void IPFPGraphEditDistance<NodeAttribute, EdgeAttribute>::
LinearSubProblem(){
  Map<MatrixXd> m_linearSubProblem(this->linearSubProblem,this->_n+1,this->_m+1);

  Map<MatrixXd> m_XkD(this->XkD,this->_n+1,this->_m+1);
  Map<MatrixXd> m_C(this->C,this->_n+1,this->_m+1);

  m_linearSubProblem = 2*m_XkD + m_C;

}



template<class NodeAttribute, class EdgeAttribute>
double IPFPGraphEditDistance<NodeAttribute, EdgeAttribute>::
getCost(double * Matrix, unsigned int n, unsigned int m){
  return linearCost(this->XkD,Matrix,n+1,m+1)+ this->Lterm;
}

template<class NodeAttribute, class EdgeAttribute>
double IPFPGraphEditDistance<NodeAttribute, EdgeAttribute>::
getCost(unsigned int * G1_to_G2,unsigned int * G2_to_G1, unsigned int n, unsigned int m){

  return linearCost(this->XkD,G1_to_G2, G2_to_G1,n,m)+ this->Lterm;
}

template<class NodeAttribute, class EdgeAttribute>
double IPFPGraphEditDistance<NodeAttribute, EdgeAttribute>::
getAlpha(){
  return this->R.back() - 2 * this->S[this->k] + this->oldLterm;
}

template<class NodeAttribute, class EdgeAttribute>
double IPFPGraphEditDistance<NodeAttribute, EdgeAttribute>::
getBeta(){
  return this->S.back() + this->S[this->k] - this->R.back() - this->oldLterm;
}


template<class NodeAttribute, class EdgeAttribute>
double * IPFPGraphEditDistance<NodeAttribute, EdgeAttribute>::
mappingsToMatrix(unsigned int * G1_to_G2,unsigned int * G2_to_G1, unsigned int n, unsigned int m, double * Matrix){
  memset(Matrix,0,sizeof(double)*(n+1)*(m+1));
  for (unsigned int i =0;i<n;i++)
    Matrix[sub2ind(i,G1_to_G2[i],n+1)] = 1;
  for (unsigned int j =0;j<m;j++)
    if (G2_to_G1[j] >= n)
      Matrix[sub2ind(G2_to_G1[j],j,n+1)] = 1;
  return Matrix;
}



template<class NodeAttribute, class EdgeAttribute>
double IPFPGraphEditDistance<NodeAttribute, EdgeAttribute>::
mappingCost( Graph<NodeAttribute,EdgeAttribute> * g1,
             Graph<NodeAttribute,EdgeAttribute> * g2,
             unsigned int* G1_to_G2 , unsigned int* G2_to_G1)
{
  return this->GedFromMapping(g1, g2, G1_to_G2, g1->Size(), G2_to_G1, g2->Size());
}

#endif // __IPFPGRAPHEDITDISTANCE_H__

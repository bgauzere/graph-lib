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

#ifndef __IPFPGRAPHEDITDISTANCE_H__
#define __IPFPGRAPHEDITDISTANCE_H__
#include <Eigen/Dense>
using namespace Eigen;
#include "hungarian-lsape.hh"
#include "GraphEditDistance.h"
#include "RandomWalksGraphEditDistance.h"
#include "utils.h"

template<class NodeAttribute, class EdgeAttribute>
class IPFPGraphEditDistance:
  public GraphEditDistance<NodeAttribute, EdgeAttribute>{
protected:
  
private:
  int maxIter = 20;
  
  GraphEditDistance<NodeAttribute,EdgeAttribute> * _ed_init;

  //Data inherent to *one* computation
  double * C;
  double * linearSubProblem;
  double * XkD;
  double * Xk;
  double Lterm;
  double oldLterm;
  double * Xkp1tD;
  double * bkp1;
  int n;
  int m;
  int k;
  std::vector<double> S;
  std::vector<double> R;
  
  void (*_MappingInit)(Graph<NodeAttribute,EdgeAttribute> * g1,
		       Graph<NodeAttribute,EdgeAttribute> * g2,
		       int * G1_to_G2, int * G2_to_G2);

  void NodeCostMatrix(Graph<NodeAttribute,EdgeAttribute> * g1,
		      Graph<NodeAttribute,EdgeAttribute> * g2);

  double * QuadraticTerm(Graph<NodeAttribute,EdgeAttribute> * g1,
			 Graph<NodeAttribute,EdgeAttribute> * g2,
			 int * G1_to_G2, int * G2_to_G1, double * XkD);
  double * QuadraticTerm(Graph<NodeAttribute,EdgeAttribute> * g1,
			 Graph<NodeAttribute,EdgeAttribute> * g2,
			 double * Matrix, double * XkD);
  double * QuadraticTerm(Graph<NodeAttribute,EdgeAttribute> * g1,
			 Graph<NodeAttribute,EdgeAttribute> * g2,
			 std::vector<std::pair<std::pair<int,int>,double> > mappings, double * XkD);

  //This linearCost is efficient for sparse Xk matrices
  double linearCost(double * CostMatrix, int * G1_to_G2,int * G2_to_G1, int n, int m);
  double linearCost(double * CostMatrix, double * Xk, int n, int m);

  // Fill this->linearSubProblem with appropriatelinear problem
  void LinearSubProblem();
  double getCost(int * G1_to_G2,int * G2_to_G1, int n, int m);
  double getAlpha();
  double getBeta();

  


public:
  IPFPGraphEditDistance(EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction,
			GraphEditDistance<NodeAttribute,EdgeAttribute> * ed_init):
    GraphEditDistance<NodeAttribute,EdgeAttribute>(costFunction),_ed_init(ed_init){};

  virtual void getOptimalMapping(Graph<NodeAttribute,EdgeAttribute> * g1,
				 Graph<NodeAttribute,EdgeAttribute> * g2,
				 int * G1_to_G2, int * G1_to_G2G2_to_G1);
  
  ~IPFPGraphEditDistance(){
    if (C != NULL)
      delete [] C;
  }

};
template<class NodeAttribute, class EdgeAttribute>
void IPFPGraphEditDistance<NodeAttribute, EdgeAttribute>::NodeCostMatrix(Graph<NodeAttribute,EdgeAttribute> * g1,
									 Graph<NodeAttribute,EdgeAttribute> * g2){
  int n=g1->Size();
  int m=g2->Size();
  
  this->C = new double[(n+1)*(m+1)];
  C[sub2ind(n,m,n)] = 0;
  for(int i=0;i<n;i++)
    for(int j=0;j<m;j++)
      C[sub2ind(i,j,(n+1))] = this->cf->NodeSubstitutionCost((*g1)[i],(*g2)[j],g1,g2);
  for(int i=0;i<n;i++)
    C[sub2ind(i,m,(n+1))] = this->cf->NodeDeletionCost((*g1)[i],g1);
  
  for(int j=0;j<m;j++)
    C[sub2ind(n,j,(n+1))] = this->cf->NodeInsertionCost((*g2)[j],g2);
}


template<class NodeAttribute, class EdgeAttribute>
double * IPFPGraphEditDistance<NodeAttribute,
			       EdgeAttribute>::QuadraticTerm(Graph<NodeAttribute,EdgeAttribute> * g1,
							     Graph<NodeAttribute,EdgeAttribute> * g2,
							     double * Matrix, double * XkD){
  int n = g1->Size();
  int m = g2->Size();
  
  std::vector<std::pair<std::pair<int,int>, double>> mappings;
  for(int i=0;i<n+1;i++)
    for(int j=0;j<m+1;j++){
      double value = Matrix[sub2ind(i,j,n+1)];
      if(value > 0.){
	std::pair<int,int> tmp = std::pair<int,int>(i,j);
	mappings.push_back(std::pair<std::pair<int,int>,double>(tmp,value));
      }
    }
  if (! XkD)
    XkD=new double[(n+1)*(m+1)];

  
  return this->QuadraticTerm(g1,g2,mappings, XkD);
    
}



template<class NodeAttribute, class EdgeAttribute>
double * IPFPGraphEditDistance<NodeAttribute,
			       EdgeAttribute>::QuadraticTerm(Graph<NodeAttribute,EdgeAttribute> * g1,
							     Graph<NodeAttribute,EdgeAttribute> * g2,
							     int * G1_to_G2, int * G2_to_G1,double * XkD){

  int n = g1->Size();
  int m = g2->Size();
  
  //Reconstruction d'un mapping
  std::vector<std::pair<std::pair<int,int>, double>> mappings;
  for (int i =0;i<n;i++){
    std::pair<int,int> tmp = std::pair<int,int>(i,G1_to_G2[i]);
    mappings.push_back(std::pair<std::pair<int,int>,double>(tmp,1.));
  }
  for (int j =0;j<m;j++)
    if (G2_to_G1[j] >= n){
      std::pair<int,int> tmp = std::pair<int,int>(G2_to_G1[j],j);
      mappings.push_back(std::pair<std::pair<int,int>,double>(tmp,1.));
    }
  if (! XkD)
    XkD=new double[(n+1)*(m+1)];
  
  return this->QuadraticTerm(g1,g2,mappings,XkD);
  
}

template<class NodeAttribute, class EdgeAttribute>
double * IPFPGraphEditDistance<NodeAttribute,
			       EdgeAttribute>::QuadraticTerm(Graph<NodeAttribute,EdgeAttribute> * g1,
							     Graph<NodeAttribute,EdgeAttribute> * g2,
							     std::vector<std::pair<std::pair<int,int>,double> > mappings,
							     double * quadraticTerm){
  int n = g1->Size();
  int m = g2->Size();
  
  if (! quadraticTerm)
    quadraticTerm=new double[(n+1)*(m+1)];

  memset(quadraticTerm,0,sizeof(double)*(n+1)*(m+1));

  for(int j = 0; j < n+1; j++){ // Attention : dans le papier sspr, condition sur x_jl /= 0. En effet, inutile pour le cas ou on multiplie a droite par le mapping. Mais nécessaire quand on utilise XtD dans le sous probleme
    for(int l = 0; l < m+1;l++){
      
      std::vector<std::pair<std::pair<int,int>,double> >::iterator it = mappings.begin();
      for(;it != mappings.end();it++){
	int i = it->first.first;
	int k = it->first.second;
	bool eps_i,eps_j,eps_k,eps_l; 
	eps_i = (i >= n);eps_j = (j >= n);eps_k = (k >= m);eps_l = (l >= m);
	   
	bool delta_e1 = false;
	if(!eps_i)
	  delta_e1 = (g1->getEdge(i,j) != NULL); // false if j>n
	bool delta_e2 = false;
	if(! eps_k)
	  delta_e2 = (g2->getEdge(k,l) != NULL);// false if l>m
	double cost = 0.0;
	//TODO : Optimize if sequence
	//If (i,j) and (k,l) are both same nodes,
	//no edges between them, so delta_e1 and delta_e2 are both 0, and so the cost
	double ced = 0.0;
	if(! eps_i)
	  ced =this->cf->EdgeDeletionCost(g1->getEdge(i,j),g1);
	double cei = 0.0;
	if(! eps_k)
	  cei = this->cf->EdgeInsertionCost(g2->getEdge(k,l),g2);

	if( ((i != j) || eps_i) && ((k != l) || eps_k)){
	  if ( (!eps_i) && (!eps_j) && (!eps_k) && (!eps_l)){
	    if (delta_e1 && delta_e2) // sub
	      cost = this->cf->EdgeSubstitutionCost(g2->getEdge(k,l), g1->getEdge(i,j),g2,g1);
	    else if ((delta_e1) && (!delta_e2)) //deletion
	      cost = ced;
	    else if ((! delta_e1) && delta_e2)
	      cost = cei;
	  }	  
	  //l is epsilon, (k,l) do not exist => deletion of (i,j)
	  else if (i < n && j < n && k < m && l >= m){
	    cost = ced*delta_e1;
	  }
	  //k is epsilon => (k,l) do not exist so delete (i,j) if exists
	  else if (i < n && j < n && k >= m && l < m){
	    cost = ced*delta_e1;
	  }
	  //(k,l) do not exists, del (i,j) if exists (factorizable with previous case)
	  else if (i < n && j < n && k >= m && l >= m){
	    cost = ced*delta_e1;
	  }
	  //i is epsilon => add (k,l) if exists
	  else if (i >= n && j < n && k < m && l < m){
	    cost = cei*delta_e2;
	  }
	  //j is epsilon => add (k,l) if exists
	  else if (i < n && j >= n && k < m && l < m){
	    cost = cei*delta_e2;
	  }
	  //i,j are epsilon, (i,j) do not exists, add (k,l) if it exists
	  else if (i >= n && j >= n && k < m && l < m){
	    cost = cei*delta_e2;
	  }
	  quadraticTerm[sub2ind(j,l,n+1)] +=  cost*it->second;
	
	}
      }
    }
  }
  return quadraticTerm;
  
 
}


template<class NodeAttribute, class EdgeAttribute>
void IPFPGraphEditDistance<NodeAttribute, EdgeAttribute>::
getOptimalMapping(Graph<NodeAttribute,EdgeAttribute> * g1,
		  Graph<NodeAttribute,EdgeAttribute> * g2,
		  int * G1_to_G2, int * G2_to_G1){
  
  this->n = g1->Size();
  this->m = g2->Size();

  S.clear();
  R.clear();
  
  NodeCostMatrix(g1,g2);
  Map<MatrixXd> m_C(this->C,n+1,m+1);
  
  //Compute Mapping init
  this->_ed_init->getOptimalMapping(g1,g2,G1_to_G2, G2_to_G1);
  
  this->XkD = this->QuadraticTerm(g1,g2,G1_to_G2, G2_to_G1,NULL);
  Map<MatrixXd> m_XkD(XkD,n+1,m+1);
  
  this->Xk = new double[(n+1)*(m+1)];
  Map<MatrixXd> m_Xk(Xk,n+1,m+1);
  
  Lterm = linearCost(this->C,G1_to_G2, G2_to_G1,n,m); 
  S.push_back(getCost(G1_to_G2, G2_to_G1,n,m));

  k=0;
  this->linearSubProblem = new double [(n+1) * (m+1)];
  Map<MatrixXd> m_linearSubProblem(linearSubProblem,n+1,m+1);

  double *u = new double[n+1];
  double *v = new double[m+1];
  //double Lterm_new = 0.0;

  this->Xkp1tD = new double [(n+1) * (m+1)];
  Map<MatrixXd> m_Xkp1tD (this->Xkp1tD,n+1,m+1);

  this->bkp1 = new double [(n+1) * (m+1)];
  Map<MatrixXd> m_bkp1 (this->bkp1,n+1,m+1);

  bool flag_continue = true;
  while((k < this->maxIter) && flag_continue){ //TODO : fixer un epsilon, param ?
    if(k!= 0)
      XkD = QuadraticTerm(g1,g2,Xk,XkD);
    else{ 
      XkD = QuadraticTerm(g1,g2,G1_to_G2,G2_to_G1,XkD);
      memset(Xk,0,sizeof(double)*(n+1)*(m+1));
      for (int i =0;i<n;i++)
	Xk[sub2ind(i,G1_to_G2[i],n+1)] = 1;
      for (int j =0;j<m;j++)
	if (G2_to_G1[j] >= n)
	  Xk[sub2ind(G2_to_G1[j],j,n+1)] = 1;
    }
    
    
    this->LinearSubProblem();//    2*m_XkD + m_C;
    hungarianLSAPE(linearSubProblem,n+1,m+1, G1_to_G2,G2_to_G1, u,v,false);

    R.push_back(linearCost(linearSubProblem,G1_to_G2, G2_to_G1,n,m));
    this->oldLterm = Lterm;
    this->Lterm= linearCost(this->C,G1_to_G2, G2_to_G1,n,m);
    XkD = QuadraticTerm(g1,g2,G1_to_G2, G2_to_G1,XkD);
    S.push_back(this->getCost(G1_to_G2, G2_to_G1,n,m));

    double alpha = getAlpha();
    double beta= getBeta();
    double t0 =0.0;
    if(beta > 0.000001)
      t0 = -alpha / (2.*beta);
    //Built a new Xk matrix (possibly not permutation)

    memset(bkp1,0,sizeof(double)*((n+1)*(m+1)));
    for(int i = 0;i<n;i++)
      bkp1[sub2ind(i,G1_to_G2[i],n+1)] = 1.;
    for(int j = 0;j<m;j++)
      if(G2_to_G1[j]>=n)      
	bkp1[sub2ind(G2_to_G1[j],j,n+1)] = 1.;
    
    if ((beta < 0.00001) || (t0 >= 1) || ( k +1 == maxIter)){
      //Check if Xk evolves
      if((m_Xk - m_bkp1).norm() < 0.0001)
	flag_continue = false;
      else
	memcpy(Xk,bkp1,sizeof(double)*(n+1)*(m+1));
      
      //Lterm = Lterm_new;
    }else{
      //Line search
      MatrixXd maj_matrix(n+1,m+1);
      maj_matrix = t0*(m_bkp1 - m_Xk);
      if(maj_matrix.norm() < 0.0001)
	flag_continue = false;
      else{
	m_Xk = m_Xk + t0*(m_bkp1 - m_Xk); 
	S[k+1] = S[k] - ((pow(alpha,2))/(4*beta));
	this->Lterm = linearCost(this->C, Xk, n,m);
      }
    }
    this->k++;    
  }
  m_Xk= m_Xk *-1;
  hungarianLSAPE(this->Xk,n+1,m+1, G1_to_G2,G2_to_G1, u,v,false);
  
  delete [] Xkp1tD;
  delete [] linearSubProblem;
  delete [] XkD;
  delete this->C;
  //  delete ed;
  delete [] u;
  delete [] v;
}



template<class NodeAttribute, class EdgeAttribute>
double IPFPGraphEditDistance<NodeAttribute, EdgeAttribute>::
linearCost(double * CostMatrix, int * G1_to_G2,int * G2_to_G1, int n, int m){
  double sum = 0.0;
  for(int i=0;i<n;i++)
    sum += CostMatrix[sub2ind(i,G1_to_G2[i],n+1)];
  for(int j=0;j<m;j++)
    if(G2_to_G1[j] >= n)
      sum+= CostMatrix[sub2ind(G2_to_G1[j],j,n+1)];
  return sum;
}

template<class NodeAttribute, class EdgeAttribute>
double IPFPGraphEditDistance<NodeAttribute, EdgeAttribute>::
linearCost(double * CostMatrix, double * X, int n, int m){
  double sum = 0.0;
  for(int i=0;i<n+1;i++)
    for(int j=0;j<m+1;j++)    
      sum += CostMatrix[sub2ind(i,j,n+1)] * X[sub2ind(i,j,n+1)];

  return sum;
}

template<class NodeAttribute, class EdgeAttribute>
void IPFPGraphEditDistance<NodeAttribute, EdgeAttribute>::
LinearSubProblem(){
  Map<MatrixXd> m_linearSubProblem(this->linearSubProblem,this->n+1,this->m+1);

  Map<MatrixXd> m_XkD(this->XkD,this->n+1,this->m+1);
  Map<MatrixXd> m_C(this->C,this->n+1,this->m+1);
  
  m_linearSubProblem = 2*m_XkD + m_C;
  
}



template<class NodeAttribute, class EdgeAttribute>
double IPFPGraphEditDistance<NodeAttribute, EdgeAttribute>::
getCost(int * G1_to_G2,int * G2_to_G1, int n, int m){
  return linearCost(this->XkD,G1_to_G2, G2_to_G1,n,m)+ this->Lterm;
}

template<class NodeAttribute, class EdgeAttribute>
double IPFPGraphEditDistance<NodeAttribute, EdgeAttribute>::
getAlpha(){
  return this->R.back() - 2 * this->S[k] + this->oldLterm;
}

template<class NodeAttribute, class EdgeAttribute>
double IPFPGraphEditDistance<NodeAttribute, EdgeAttribute>::
getBeta(){
  return S.back() + S[k] -R.back() - this->oldLterm;
}



#endif // __IPFPGRAPHEDITDISTANCE_H__

/**
 * @file RiesenCostMatrix.h
 * @author bgauzere <<benoit.gauzere@insa-rouen.fr>> 
 * @version     0.0.1 - Fri Apr 13 2018
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 */

#ifndef __RIESENCOSTMATRIX_H__
#define __RIESENCOSTMATRIX_H__

#include "CostMatrix.h"
#include "gl_utils.h"
#include <lsape.h>
template<class NodeAttribute, class EdgeAttribute>
class RiesenCostMatrix:
  public CostMatrix<NodeAttribute, EdgeAttribute>
{
protected:
  
  enum lsape::LSAPE_MODEL  my_solver; 
  double SubstitutionCost(GNode<NodeAttribute,EdgeAttribute> * v1,
			  GNode<NodeAttribute,EdgeAttribute> * v2,
			  Graph<NodeAttribute,EdgeAttribute> * g1,
			  Graph<NodeAttribute,EdgeAttribute> * g2);
  
  double DeletionCost(GNode<NodeAttribute,EdgeAttribute> * v1,
		      Graph<NodeAttribute,EdgeAttribute> * g1);

  double InsertionCost(GNode<NodeAttribute,EdgeAttribute> * v2,
		       Graph<NodeAttribute,EdgeAttribute> * g2);
public:
  RiesenCostMatrix(EditDistanceCost<NodeAttribute,EdgeAttribute> * cf):
    CostMatrix<NodeAttribute, EdgeAttribute>(cf),my_solver(lsape::ECBP){};

  void computeCostMatrix(Graph<NodeAttribute,EdgeAttribute> * g1,
			 Graph<NodeAttribute,EdgeAttribute> * g2);
  
  
};

template<class NodeAttribute, class EdgeAttribute>
void RiesenCostMatrix<NodeAttribute, EdgeAttribute>::
computeCostMatrix(Graph<NodeAttribute,EdgeAttribute> * g1,
		  Graph<NodeAttribute,EdgeAttribute> * g2){

  if(this->_C)
    delete [] this->_C;

  int n=g1->Size();
  int m=g2->Size();
  this->_n = n;
  this->_m = m;
  this->_C = new double[(n+1) * (m+1)];
  this->_C = (double*)memset(this->_C,0,sizeof(double)*(n+1) * (m+1));
  for (int i =0;i<n;i++)
    for(int j = 0;j<m;j++)
      this->_C[sub2ind(i,j,n+1)] = this->SubstitutionCost((*g1)[i],(*g2)[j],g1,g2);
  
  for (int i =0;i<n;i++)
    this->_C[sub2ind(i,m,n+1)] = this->DeletionCost((*g1)[i],g1);

  for (int j =0;j<m;j++)
    this->_C[sub2ind(n,j,n+1)] = this->InsertionCost((*g2)[j],g2);

  this->_C[sub2ind(n,m,n+1)] = 0;
};
  
template<class NodeAttribute, class EdgeAttribute>
double RiesenCostMatrix<NodeAttribute, EdgeAttribute>::
SubstitutionCost(GNode<NodeAttribute,EdgeAttribute> * v1,
		 GNode<NodeAttribute,EdgeAttribute> * v2,
		 Graph<NodeAttribute,EdgeAttribute> * g1,
		 Graph<NodeAttribute,EdgeAttribute> * g2){
  int n=v1->Degree();
  int m=v2->Degree();

  GEdge<EdgeAttribute> * e1 = v1->getIncidentEdges(); //edge from v1 in G1
  GEdge<EdgeAttribute> * _e2 = v2->getIncidentEdges(); //edge from v2 in G2

  GEdge<EdgeAttribute> * e2 = _e2; //We keep a copy of e2 to start again an iteration

  double * local_C = new double[(n+1) * (m+1)];
  memset(local_C,0,sizeof(double)*(n+1) * (m+1));
  for (int i=0;e1;i++){
    e2 = _e2; 
    for (int j=0;e2;j++){
      local_C[sub2ind(i,j,n+1)] = this->_cf->EdgeSubstitutionCost(e1,e2,g1,g2);
      // this->cf->NodeSubstitutionCost((*g1)[e1->IncidentNode()],
      // 								 (*g2)[e2->IncidentNode()],
      // 								 g1,g2) +
      e2 = e2->Next();
    }
    e1 = e1->Next();
  }

  e1 = v1->getIncidentEdges();
  for (int i=0;e1;i++){
    local_C[sub2ind(i,m,n+1)] = this->_cf->EdgeDeletionCost(e1,g1); //this->cf->NodeDeletionCost((*g1)[e1->IncidentNode()],g1) +
    e1 = e1->Next();
  }
  e2 = v2->getIncidentEdges();
  for (int j=0;e2;j++){
    local_C[sub2ind(n, j, n+1)] = this->_cf->EdgeInsertionCost(e2,g2); //this->cf->NodeInsertionCost((*g2)[e2->IncidentNode()],g2) +
    e2 = e2->Next();
  }
  local_C[sub2ind(n,m,n+1)] = 0;
  unsigned int *rho = new unsigned int[n+1];
  unsigned int *varrho = new unsigned int[m+1];
  double *u = new double[n+1];
  double *v = new double[m+1];

  rho = (unsigned int*)memset((void*)rho,0,sizeof(unsigned int)*(n+1));
  varrho = (unsigned int*)memset((void*)varrho,0,sizeof(unsigned int)*(m+1));
  u = (double*)memset((void*)u,0,sizeof(double)*(n+1));
  v = (double*)memset((void*)v,0,sizeof(double)*(m+1));

  lsape::lsapeSolver<double>(local_C,n+1,m+1, rho, varrho, u,v,this->my_solver,1);
  double cost=0.0;
  for (int i =0;i<n+1;i++)
    cost += u[i];
  for (int j =0;j<m+1;j++)
    cost += v[j];
  delete [] u;delete [] v;
  delete [] rho;delete [] varrho;
  delete [] local_C;
  return cost + this->_cf->NodeSubstitutionCost(v1,v2,g1,g2);
}


template<class NodeAttribute, class EdgeAttribute>
double RiesenCostMatrix<NodeAttribute, EdgeAttribute>::
DeletionCost(GNode<NodeAttribute,EdgeAttribute> * v1,Graph<NodeAttribute,EdgeAttribute> * g1){
  int n=v1->Degree();

  GEdge<EdgeAttribute> * e1 = v1->getIncidentEdges();
  double cost = 0.0;
  for (int i=0;i<n;i++, e1 = e1->Next())
    cost += this->_cf->EdgeDeletionCost(e1,g1); //this->_cf->NodeDeletionCost((*g1)[e1->IncidentNode()],g1) +
  
  return cost + this->_cf->NodeDeletionCost(v1,g1);
}

template<class NodeAttribute, class EdgeAttribute>
double RiesenCostMatrix<NodeAttribute, EdgeAttribute>::
InsertionCost(GNode<NodeAttribute,EdgeAttribute> * v2,Graph<NodeAttribute,EdgeAttribute> * g2){
  int n=v2->Degree();

  GEdge<EdgeAttribute> * e2 = v2->getIncidentEdges();
  double cost = 0.0;
  for (int i=0;i<n;i++, e2 = e2->Next())
    cost  += this->_cf->EdgeInsertionCost(e2,g2);// this->_cf->NodeInsertionCost((*g2)[e2->IncidentNode()],g2) +
  return cost + this->_cf->NodeInsertionCost(v2,g2);
}



#endif // __RIESENCOSTMATRIX_H__

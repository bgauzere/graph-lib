/**
 * @file BipartiteGraphEditDistance.h
 * @author Benoit <<benoit.gauzere@insa-rouen.fr>> 
 * @version     0.0.1 - Mon Feb  6 2017
 *
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *
 * Bipartite Graph edit distance algorithm as described in Riesen's
 * book.  Riesen, K. (2015). Structural pattern recognition with graph
 * edit distance. Advances in Computer Vision and Pattern Recognition,
 * Cham.
 */

#ifndef __BIPARTITEGRAPHEDITDISTANCE_H__
#define __BIPARTITEGRAPHEDITDISTANCE_H__

#include <lsape.h>
#include "GraphEditDistance.h"
#include "gl_utils.h"



//TODO : donner la possibilité de récupérer le mapping ?
template<class NodeAttribute, class EdgeAttribute>
class BipartiteGraphEditDistance:
  public GraphEditDistance<NodeAttribute, EdgeAttribute>
{
private:
  enum lsape::LSAPE_MODEL  my_solver;

protected:
  double * C;

protected:
  virtual void computeCostMatrix(Graph<NodeAttribute,EdgeAttribute> * g1,
				 Graph<NodeAttribute,EdgeAttribute> * g2);
         
  double SubstitutionCost(GNode<NodeAttribute,EdgeAttribute> * v1,
			  GNode<NodeAttribute,EdgeAttribute> * v2,
			  Graph<NodeAttribute,EdgeAttribute> * g1,
			  Graph<NodeAttribute,EdgeAttribute> * g2);

  double DeletionCost(GNode<NodeAttribute,EdgeAttribute> * v1,
		      Graph<NodeAttribute,EdgeAttribute> * g1);

  double InsertionCost(GNode<NodeAttribute,EdgeAttribute> * v2,
		       Graph<NodeAttribute,EdgeAttribute> * g2);

public:
  BipartiteGraphEditDistance(EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction):
    GraphEditDistance<NodeAttribute,EdgeAttribute>(costFunction),my_solver(lsape::ECBP),
    C(NULL)
  {};

  BipartiteGraphEditDistance(EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction,
			     enum lsape::LSAPE_MODEL f):
    GraphEditDistance<NodeAttribute,EdgeAttribute>(costFunction),
    my_solver(f),C(NULL){};

  void setSolver(enum lsape::LSAPE_MODEL f){
    this->my_solver = f;
  }

  // virtual double operator()(Graph<NodeAttribute,EdgeAttribute> * g1,
  // 			    Graph<NodeAttribute,EdgeAttribute> * g2);
  double getTimedOptimalMapping(Graph<NodeAttribute,EdgeAttribute> * g1,
				Graph<NodeAttribute,EdgeAttribute> * g2,
				unsigned int * G1_to_G2, unsigned int * G2_to_G1);
  virtual void getOptimalMapping(Graph<NodeAttribute,EdgeAttribute> * g1,
				 Graph<NodeAttribute,EdgeAttribute> * g2,
				 unsigned int * G1_to_G2, unsigned int * G2_to_G2);
  double getUpperBound(Graph<NodeAttribute,EdgeAttribute> * g1,
		       Graph<NodeAttribute,EdgeAttribute> * g2,
		       double & time);
  

    virtual ~BipartiteGraphEditDistance(){
    if (this->C != NULL) delete [] this->C;
  }

  virtual BipartiteGraphEditDistance<NodeAttribute, EdgeAttribute> * clone() const
  {
    return new BipartiteGraphEditDistance<NodeAttribute, EdgeAttribute> (*this);
  }

};


//A remonter a GraphEditDistance.h ?
template<class NodeAttribute, class EdgeAttribute>
double BipartiteGraphEditDistance<NodeAttribute, EdgeAttribute>::
getTimedOptimalMapping(Graph<NodeAttribute,EdgeAttribute> * g1,
		       Graph<NodeAttribute,EdgeAttribute> * g2,
		       unsigned int * G1_to_G2, unsigned int * G2_to_G1){
  int n=g1->Size();
  int m=g2->Size();
  // Compute C
  delete [] this->C;
  computeCostMatrix(g1,g2);
  //Compute optimal assignement
  double *u = new double[n+1];
  double *v = new double[m+1];
  clock_t t = clock();
  lsape::lsapeSolver<double>(C,n+1,m+1, G1_to_G2, G2_to_G1, u,v,my_solver,1);
  t = clock() - t;
  delete [] u;
  delete [] v;
  return ((float)t) / CLOCKS_PER_SEC;

}

template<class NodeAttribute, class EdgeAttribute>
void BipartiteGraphEditDistance<NodeAttribute, EdgeAttribute>::
getOptimalMapping(Graph<NodeAttribute,EdgeAttribute> * g1,
		  Graph<NodeAttribute,EdgeAttribute> * g2,
		  unsigned int * G1_to_G2, unsigned int * G2_to_G1){
  int n=g1->Size();
  int m=g2->Size();
  // Compute C
  delete [] this->C;
  computeCostMatrix(g1,g2);
  //Compute optimal assignement
  double *u = new double[n+1];
  double *v = new double[m+1];
  lsape::lsapeSolver<double>(C,n+1,m+1, G1_to_G2, G2_to_G1, u,v,my_solver,1);
  delete [] u;
  delete [] v;

}

template<class NodeAttribute, class EdgeAttribute>
double BipartiteGraphEditDistance<NodeAttribute, EdgeAttribute>::
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
      local_C[sub2ind(i,j,n+1)] = this->cf->EdgeSubstitutionCost(e1,e2,g1,g2);
      // this->cf->NodeSubstitutionCost((*g1)[e1->IncidentNode()],
      // 								 (*g2)[e2->IncidentNode()],
      // 								 g1,g2) +
      e2 = e2->Next();
    }
    e1 = e1->Next();
  }

  e1 = v1->getIncidentEdges();
  for (int i=0;e1;i++){
    local_C[sub2ind(i,m,n+1)] = this->cf->EdgeDeletionCost(e1,g1); //this->cf->NodeDeletionCost((*g1)[e1->IncidentNode()],g1) +
    e1 = e1->Next();
  }
  e2 = v2->getIncidentEdges();
  for (int j=0;e2;j++){
    local_C[sub2ind(n, j, n+1)] = this->cf->EdgeInsertionCost(e2,g2); //this->cf->NodeInsertionCost((*g2)[e2->IncidentNode()],g2) +
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

  lsape::lsapeSolver<double>(local_C,n+1,m+1, rho, varrho, u,v,my_solver,1);
  double cost=0.0;
  for (int i =0;i<n+1;i++)
    cost += u[i];
  for (int j =0;j<m+1;j++)
    cost += v[j];
  delete [] u;delete [] v;
  delete [] rho;delete [] varrho;
  delete [] local_C;
  return cost + this->cf->NodeSubstitutionCost(v1,v2,g1,g2);
}


template<class NodeAttribute, class EdgeAttribute>
double BipartiteGraphEditDistance<NodeAttribute, EdgeAttribute>::
DeletionCost(GNode<NodeAttribute,EdgeAttribute> * v1,Graph<NodeAttribute,EdgeAttribute> * g1){
  int n=v1->Degree();

  GEdge<EdgeAttribute> * e1 = v1->getIncidentEdges();
  double cost = 0.0;
  for (int i=0;i<n;i++, e1 = e1->Next())
    cost += this->cf->EdgeDeletionCost(e1,g1); //this->cf->NodeDeletionCost((*g1)[e1->IncidentNode()],g1) +
  
  return cost + this->cf->NodeDeletionCost(v1,g1);
}

template<class NodeAttribute, class EdgeAttribute>
double BipartiteGraphEditDistance<NodeAttribute, EdgeAttribute>::
InsertionCost(GNode<NodeAttribute,EdgeAttribute> * v2,Graph<NodeAttribute,EdgeAttribute> * g2){
  int n=v2->Degree();

  GEdge<EdgeAttribute> * e2 = v2->getIncidentEdges();
  double cost = 0.0;
  for (int i=0;i<n;i++, e2 = e2->Next())
    cost  += this->cf->EdgeInsertionCost(e2,g2);// this->cf->NodeInsertionCost((*g2)[e2->IncidentNode()],g2) +
  return cost + this->cf->NodeInsertionCost(v2,g2);
}

template<class NodeAttribute, class EdgeAttribute>
void BipartiteGraphEditDistance<NodeAttribute, EdgeAttribute>::
computeCostMatrix(Graph<NodeAttribute,EdgeAttribute> * g1,
		  Graph<NodeAttribute,EdgeAttribute> * g2){

  //delete [] C;

  int n=g1->Size();
  int m=g2->Size();
  // if (C != NULL){
  //   delete [] C;
  //   C = NULL;
  // }
  C = new double[(n+1) * (m+1)];
  C = (double*)memset(C,0,sizeof(double)*(n+1) * (m+1));
  for (int i =0;i<n;i++)
    for(int j = 0;j<m;j++)
      C[sub2ind(i,j,n+1)] = this->SubstitutionCost((*g1)[i],(*g2)[j],g1,g2);
  
  for (int i =0;i<n;i++)
    C[sub2ind(i,m,n+1)] = this->DeletionCost((*g1)[i],g1);

  for (int j =0;j<m;j++)
    C[sub2ind(n,j,n+1)] = this->InsertionCost((*g2)[j],g2);

  C[sub2ind(n,m,n+1)] = 0;
}

template<class NodeAttribute, class EdgeAttribute>
double BipartiteGraphEditDistance<NodeAttribute, EdgeAttribute>::
getUpperBound(Graph<NodeAttribute,EdgeAttribute> * g1,
	      Graph<NodeAttribute,EdgeAttribute> * g2,double & comptime){
  int n=g1->Size();
  int m=g2->Size();
  unsigned int * G1_to_G2 = new unsigned int[n];
  unsigned int * G2_to_G1 = new unsigned int[m];
  comptime = this->getTimedOptimalMapping(g1,g2,G1_to_G2,G2_to_G1);
  double ged = this->GedFromMapping(g1,g2,G1_to_G2,n,G2_to_G1,m);
  delete [] G1_to_G2;
  delete [] G2_to_G1;
  return ged;
}


#endif // __BIPARTITEGRAPHEDITDISTANCE_H__

/**
 * @file BipartiteLowerBoundGraphEditDistance.h
 * @author Benoit <<benoit.gauzere@insa-rouen.fr>> 
 * @version     0.0.1 - Mon Feb  6 2017
 *
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *
 * Bipartite Graph edit distance algorithm configured to compute a
 * lower bound of ged, as described in Riesen's book.  Riesen,
 * K. (2015). Structural pattern recognition with graph edit
 * distance. Advances in Computer Vision and Pattern Recognition,
 * Cham.
 */

#ifndef __BIPARTITELOWERBOUNDGRAPHEDITDISTANCE_H__
#define __BIPARTITELOWERBOUNDGRAPHEDITDISTANCE_H__

#include <lsape.h>
#include "GraphEditDistance.h"

#include "gl_utils.h"

//TODO : donner la possibilité de récupérer le mapping ?
template<class NodeAttribute, class EdgeAttribute>
class BipartiteLowerBoundEditDistance:
  public GraphEditDistance<NodeAttribute, EdgeAttribute>
{

private:
  lsape::LSAPE_MODEL my_solver;
  
protected:
  double * C;

protected:
  //Compute a n+1 \times m+1 cost matrix
  void computeCostMatrix(Graph<NodeAttribute,EdgeAttribute> * g1,
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
  BipartiteLowerBoundEditDistance(EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction):
    GraphEditDistance<NodeAttribute,EdgeAttribute>(costFunction),
    my_solver(lsape::ECBP),C(NULL){
    
  };

  BipartiteLowerBoundEditDistance(EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction,
				  lsape::LSAPE_MODEL f):
    GraphEditDistance<NodeAttribute,EdgeAttribute>(costFunction),
    my_solver(f),C(NULL){};

  void setSolver(lsape::LSAPE_MODEL f){
    this->my_solver = f;
  }


  double getTimedOptimalMapping(Graph<NodeAttribute,EdgeAttribute> * g1,
				Graph<NodeAttribute,EdgeAttribute> * g2,
				unsigned int * G1_to_G2, unsigned int * G2_to_G1, double &opt_cost);
  void getOptimalMapping(Graph<NodeAttribute,EdgeAttribute> * g1,
			 Graph<NodeAttribute,EdgeAttribute> * g2,
			 unsigned int * G1_to_G2, unsigned int * G2_to_G1);

  double getLowerBound(Graph<NodeAttribute,EdgeAttribute> * g1,
		       Graph<NodeAttribute,EdgeAttribute> * g2);

  double getLowerBound(Graph<NodeAttribute,EdgeAttribute> * g1,
		       Graph<NodeAttribute,EdgeAttribute> * g2,
		       double & time);

  ~BipartiteLowerBoundEditDistance(){
    if (this->C != NULL) delete [] this->C;
  }
  
};


//A remonter a GraphEditDistance.h ?




template<class NodeAttribute, class EdgeAttribute>
void BipartiteLowerBoundEditDistance<NodeAttribute, EdgeAttribute>::
getOptimalMapping(Graph<NodeAttribute,EdgeAttribute> * g1,
		  Graph<NodeAttribute,EdgeAttribute> * g2,
		  unsigned int * G1_to_G2,unsigned int * G2_to_G1){
  double cost;
  __attribute__((unused)) double tmp = getTimedOptimalMapping(g1,g2,G1_to_G2,G2_to_G1,cost);

}

template<class NodeAttribute, class EdgeAttribute>
double BipartiteLowerBoundEditDistance<NodeAttribute, EdgeAttribute>::
getTimedOptimalMapping(Graph<NodeAttribute,EdgeAttribute> * g1,
		       Graph<NodeAttribute,EdgeAttribute> * g2,
		       unsigned int * G1_to_G2,unsigned int * G2_to_G1, double & opt_cost){
  int n=g1->Size();
  int m=g2->Size();
  // Compute C
  delete [] this->C;
  computeCostMatrix(g1,g2);
  
  //Compute optimal assignement
  double *u = new double[n+1];
  memset(u,0,sizeof(double)*(n+1));
  double *v = new double[m+1];
  memset(v,0,sizeof(double)*(m+1));
    
  // for (int i = 0;i<n+1;i++){
  //   for (int j =0;j<m+1;j++)
  //     std::cout << C[sub2ind(i,j,n+1)] <<",";
  //   std::cout << std::endl;
  // }
  // std::cout << "Diag :" << std::endl;
  // for (int i = 0;i<n+1;i++){
  //   std::cout << C[sub2ind(i,i,n+1)] <<",";
  // }

  clock_t t = clock();
  lsape::lsapeSolver<double>(C,n+1,m+1, G1_to_G2, G2_to_G1, u,v,my_solver,1);
  t = clock() - t;
  opt_cost =0.0;
  for (int i =0;i<n+1;i++)
    opt_cost += u[i];
  for (int j =0;j<m+1;j++)
    opt_cost += v[j];

  // for (int i =0;i<n+1;i++)
  //   std::cout << u[i] << std::endl;
  // for (int j=0;j<m+1;j++)
  //   std::cout << v[j] << std::endl;

  // std::cout << opt_cost << std::endl;

  delete [] u;
  delete [] v;
  return ((float)t) / CLOCKS_PER_SEC;

    
  
}

template<class NodeAttribute, class EdgeAttribute>
double BipartiteLowerBoundEditDistance<NodeAttribute, EdgeAttribute>::
getLowerBound(Graph<NodeAttribute,EdgeAttribute> * g1,
	      Graph<NodeAttribute,EdgeAttribute> * g2, double & time){
  int n = g1->Size();
  int m = g2->Size();
  unsigned int * G1_to_G2 = new unsigned int[n+1];
  unsigned int * G2_to_G1 = new unsigned int[m+1];
  double lower_bound = 0.0;
  
  time = this->getTimedOptimalMapping(g1,g2,G1_to_G2,G2_to_G1,lower_bound);
			  
  // double lower_bound = 0.0;
  // for(int i =0;i<n;i++)
  //   lower_bound += C[sub2ind(i,G1_to_G2[i],n+1)];
  // for(int j =0;j<m;j++)
  //   if(G2_to_G1[j] == n)
  //     lower_bound += C[sub2ind(n,j,n+1)];

  // std::cout << lower_bound << std::endl;
  
  delete [] G1_to_G2;
  delete [] G2_to_G1;
  
  return lower_bound;
}
  
template<class NodeAttribute, class EdgeAttribute>
double BipartiteLowerBoundEditDistance<NodeAttribute, EdgeAttribute>::
getLowerBound(Graph<NodeAttribute,EdgeAttribute> * g1,
	      Graph<NodeAttribute,EdgeAttribute> * g2){
  __attribute__((unused)) double tmp;
  return getLowerBound(g1,g2,tmp);
}

template<class NodeAttribute, class EdgeAttribute>
double BipartiteLowerBoundEditDistance<NodeAttribute, EdgeAttribute>::
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
  rho = (unsigned int*)memset((void*)rho,0,sizeof(unsigned int)*(n+1));
  unsigned int *varrho = new unsigned int[m+1];
  varrho = (unsigned int*)memset((void*)varrho,0,sizeof(unsigned int)*(m+1));
  double *u = new double[n+1];
  u = (double*)memset((void*)u,0,sizeof(double)*(n+1));
  double *v = new double[m+1];
  v = (double*)memset((void*)v,0,sizeof(double)*(m+1));
  
  
  lsape::lsapeSolver<double>(local_C,n+1,m+1, rho, varrho, u,v,my_solver,1);
  
  double cost=0.0;
  for (int i =0;i<n+1;i++)
    cost += u[i];
  for (int j =0;j<m+1;j++)
    cost += v[j];
  delete [] u;delete [] v;
  delete [] rho;delete [] varrho;
  //delete [] local_C;
  // std::cout << cost << std::endl;
  // std::cout << this->cf->NodeSubstitutionCost(v1,v2,g1,g2) << std::endl;
    
  return 0.5*cost + this->cf->NodeSubstitutionCost(v1,v2,g1,g2);
}


template<class NodeAttribute, class EdgeAttribute>
double BipartiteLowerBoundEditDistance<NodeAttribute, EdgeAttribute>::
DeletionCost(GNode<NodeAttribute,EdgeAttribute> * v1,Graph<NodeAttribute,EdgeAttribute> * g1){
  int n=v1->Degree();

  GEdge<EdgeAttribute> * e1 = v1->getIncidentEdges();
  double cost = 0.0;
  for (int i=0;i<n;i++, e1 = e1->Next())
    cost += this->cf->EdgeDeletionCost(e1,g1); //this->cf->NodeDeletionCost((*g1)[e1->IncidentNode()],g1) +
  
  return 0.5*cost + this->cf->NodeDeletionCost(v1,g1);
}

template<class NodeAttribute, class EdgeAttribute>
double BipartiteLowerBoundEditDistance<NodeAttribute, EdgeAttribute>::
InsertionCost(GNode<NodeAttribute,EdgeAttribute> * v2,Graph<NodeAttribute,EdgeAttribute> * g2){
  int n=v2->Degree();

  GEdge<EdgeAttribute> * e2 = v2->getIncidentEdges();
  double cost = 0.0;
  for (int i=0;i<n;i++, e2 = e2->Next())
    cost  += this->cf->EdgeInsertionCost(e2,g2);// this->cf->NodeInsertionCost((*g2)[e2->IncidentNode()],g2) +
  return 0.5*cost + this->cf->NodeInsertionCost(v2,g2);
}

template<class NodeAttribute, class EdgeAttribute>
void BipartiteLowerBoundEditDistance<NodeAttribute, EdgeAttribute>::
computeCostMatrix(Graph<NodeAttribute,EdgeAttribute> * g1,
		  Graph<NodeAttribute,EdgeAttribute> * g2){
  int n=g1->Size();
  int m=g2->Size();
  this->C = new double[(n+1) * (m+1)];
  memset(this->C,0,sizeof(double)*(n+1) * (m+1));
  for (int i =0;i<n;i++)
    for(int j = 0;j<m;j++){
      // std::cout << i << "," << j << std::endl;
      this->C[sub2ind(i,j,n+1)] = this->SubstitutionCost((*g1)[i],(*g2)[j],g1,g2);
    }
  for (int i =0;i<n;i++)
    this->C[sub2ind(i,m,n+1)] = this->DeletionCost((*g1)[i],g1);

  for (int j =0;j<m;j++)
    this->C[sub2ind(n,j,n+1)] = this->InsertionCost((*g2)[j],g2);

  this->C[sub2ind(n,m,n+1)] = 0;
}
#endif // __BIPARTITEGRAPHEDITDISTANCE_H__

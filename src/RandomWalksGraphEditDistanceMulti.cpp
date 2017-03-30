/*
 * @file RandomWalksGraphEditDistanceMulti.cpp
 * @author Evariste <<evariste.daller@unicaen.fr>>
 * @version     0.0.1 - Wed Mar 29 2017
 */

 #include "RandomWalksGraphEditDistanceMulti.h"
 #include "AllPerfectMatchings-ec.h"

std::list<int*> RandomWalksGraphEditDistanceMulti::
getKOptimalMappings(Graph<int,int> * g1,
                    Graph<int,int> * g2,
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



 double RandomWalksGraphEditDistanceMulti::
 operator() (Graph<int,int> * g1,
             Graph<int,int> * g2,
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
    std::cout << nged << std::endl;

    if (ged > nged || ged == -1)
      ged = nged;
  }

  for (it=mappings.begin(); it!=mappings.end(); it++)
    delete[] *it;

  std::cout << std::endl;
  delete [] G2_to_G1;
  return ged;
}




 double RandomWalksGraphEditDistanceMulti::
 operator() (Graph<int,int> * g1,
             Graph<int,int> * g2)
{
  return operator() (g1,g2,nep);
}

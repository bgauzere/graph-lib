/*
 * @file RandomWalksGraphEditDistanceMulti.cpp
 * @author Evariste <<evariste.daller@unicaen.fr>>
 * @version     0.0.1 - Wed Mar 29 2017
 */

 #include "RandomWalksGraphEditDistanceMulti.h"
 #include "AllPerfectMatchings-ec.h"
 #include "hungarian-lsap.hh"


 void RandomWalksGraphEditDistanceMulti::
 computeCostMatrixLSAP(int n, int m)
 {
   Clsap = new  double[(n+m)*(n+m)];
   //memset(Clsap, -1.0, sizeof(double)*(n+m)*(n+m)); // inf costs to all non-possible mappings
   for (int j=0; j<m+n; j++)
     for (int i=0; i<m+n; i++)
       if (i>=n && j>=m) Clsap[sub2ind(i,j,n+m)] = 0;
       else Clsap[sub2ind(i,j,n+m)] = -1.0;

   //XXX changer column first
   for(int i=0;i<n;i++)
     for(int j=0;j<m;j++)
       Clsap[sub2ind(i,j,n+m)] = this->C[sub2ind(i,j,n+1)];

    for(int i=0;i<n;i++)
      Clsap[sub2ind(i,m+i,n+m)] = this->C[sub2ind(i,m,n+1)];
    for(int j=0;j<m;j++)
      Clsap[sub2ind(n+j,j,n+m)] = this->C[sub2ind(n,j,n+1)];
 }


std::list<int*> RandomWalksGraphEditDistanceMulti::
getKOptimalMappings(Graph<int,int> * g1,
                    Graph<int,int> * g2,
                    const int& k)
{
  int n=g1->Size();
  int m=g2->Size();

  // Compute C
  delete [] this->C;     this->C = NULL;
  delete [] this->Clsap; this->Clsap = NULL;

#if XP_OUTPUT
  auto start = std::chrono::steady_clock::now();
#endif
  this->computeCostMatrix(g1, g2);
  this->computeCostMatrixLSAP(n,m);
#if XP_OUTPUT
  auto end = std::chrono::steady_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  std::cout << elapsed.count() << ":";

  start = std::chrono::steady_clock::now();
#endif

// the returned mappings
std::list<int*> mappings;

#if XP_OUTPUT
  // perform XP_TIME_SAMPLES mesures
  for (int _nts=0; _nts<XP_TIME_SAMPLES; _nts++){
  mappings.clear();
#endif

  // Compute an optimal assignement
  double *u = new double[n+1];
  double *v = new double[m+1];
  int* G1_to_G2 = new int[n+1];
  int* G2_to_G1 = new int[m+1];

  //hungarianLSAP<double,int>(this->Clsap, n+m, m+n, G1_to_G2, u, v, true);
  hungarianLSAPE(this->C, n+1, m+1, G1_to_G2, G2_to_G1, u, v, false);

  // Compute LSAP solution from LSAPE
  int* rhoperm = new int[n+m];
  bool* epsAssign = new bool[n]; // is eps[i] assigned
  for (int j=0; j<n; j++) epsAssign[j] = false;
  for (int i=0; i<n; i++){
    if (G1_to_G2[i] < m)
      rhoperm[i] = G1_to_G2[i];
    else{
      rhoperm[i] = i+m;
      epsAssign[i] = true;
    }
  }
  int firstEpsNonAssign = 0;
  for (int j=0; j<m; j++){
    if (G2_to_G1[j] == n)
      rhoperm[j+n] = j;
    //*
    else{ // find the first epsilon not assigned
      while (firstEpsNonAssign < n && epsAssign[firstEpsNonAssign]) firstEpsNonAssign++;
      rhoperm[j+n] = firstEpsNonAssign + m;
      epsAssign[firstEpsNonAssign] = true;
    }
    //*/
  }

  // Compute LSAP u and v
  double *lu = new double[n+m];
  double *lv = new double[n+m];
  for (int i=0; i<n; i++) lu[i] = u[i];
  for (int i=n; i<n+m; i++) lu[i] = 0;
  for (int j=0; j<m; j++) lv[j] = v[j];
  for (int j=m; j<n+m; j++) lv[j] = 0;


  // Compute the k optimal mappings
  cDigraph<int> edg = equalityDigraph<double,int> (this->Clsap, n+m, n+m, rhoperm, lu, lv);
  AllPerfectMatchings<int> apm(edg);
  apm.enumPerfectMatchings(edg,k);
  mappings = apm.getPerfectMatchings();

  // Add the first one to the list
  mappings.push_front(rhoperm);

  delete [] epsAssign;
  delete [] u;
  delete [] v;
  delete [] lu;
  delete [] lv;
  delete [] G2_to_G1;
  delete [] G1_to_G2;

#if XP_OUTPUT
  if (_nts < XP_TIME_SAMPLES-1) delete [] rhoperm;
  } // end for 1..XP_TIME_SAMPLES
  end = std::chrono::steady_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  std::cout << elapsed.count() / XP_TIME_SAMPLES << ":";
#endif


  return mappings;
}


void RandomWalksGraphEditDistanceMulti::
getOptimalMapping (Graph<int,int> * g1,
                   Graph<int,int> * g2,
                   int * G1_to_G2, int * G2_to_G1 )
{
  int n=g1->Size();
  int m=g2->Size();

  std::list<int*> mappings = getKOptimalMappings(g1, g2, nep);

#if XP_OUTPUT
  auto start = std::chrono::steady_clock::now();
#endif

  // Get the min of ged;
  ged = -1.0;
  double nged;

  typename std::list<int*>::const_iterator it;
  for (it=mappings.begin(); it!=mappings.end(); it++){
    int* lsapMapping = *it;

    // computation of G1_to_G2 and G2_to_G1
    for (int j=0; j<m; j++){ // connect all to epsilon by default
      G2_to_G1[j] = n;
    }

    for (int i=0; i<n; i++){
      if (lsapMapping[i] >= m)
        G1_to_G2[i] = m; // i -> epsilon
      else{
        G1_to_G2[i] = lsapMapping[i];
        G2_to_G1[lsapMapping[i]] = i;
      }
    }

    for (int j=0; j<m; j++){
      if (lsapMapping[n+j] < m){
        G2_to_G1[j] = n; // epsilon -> j
      }
    }

    nged = this->GedFromMapping(g1, g2, G1_to_G2,n, G2_to_G1,m);

    if (ged > nged || ged == -1)
      ged = nged;
  }

#if XP_OUTPUT
  auto end = std::chrono::steady_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  std::cout << elapsed.count() << ":";
  std::cout << (mappings.size() == nep) << ":";
#endif

  for (it=mappings.begin(); it!=mappings.end(); it++)
    delete[] *it;

}


 double RandomWalksGraphEditDistanceMulti::
 operator() (Graph<int,int> * g1,
             Graph<int,int> * g2,
             const int& k)
{
  int n=g1->Size();
  int m=g2->Size();
  std::list<int*> mappings = getKOptimalMappings(g1, g2, k);

#if XP_OUTPUT
  auto start = std::chrono::steady_clock::now();
#endif

  ged = -1;  double nged;
  int* G1_to_G2 = new int[n+1];
  int* G2_to_G1 = new int[m+1];

  //XXX There is a bug here :
  //    if getKOptimalMappings is called after the allocation of G1_to_G2 and G2_to_G1
  //    then the cost matrix computed is (sometimes) false
  //    This happens even if getKOptimalMappings is called from another member function as getOptimalMapping
  //TODO fix this bug so we can avoid code pasting like this
  //this->nep=k;
  //this->getOptimalMapping(g1, g2, G1_to_G2, G2_to_G1);

  typename std::list<int*>::const_iterator it;
  for (it=mappings.begin(); it!=mappings.end(); it++){
    int* lsapMapping = *it;

    // computation of G1_to_G2 and G2_to_G1
    for (int j=0; j<m; j++){ // connect all to epsilon by default
      G2_to_G1[j] = n;
    }

    for (int i=0; i<n; i++){
      if (lsapMapping[i] >= m)
        G1_to_G2[i] = m; // i -> epsilon
      else{
        G1_to_G2[i] = lsapMapping[i];
        G2_to_G1[lsapMapping[i]] = i;
      }
    }

    for (int j=0; j<m; j++){
      if (lsapMapping[n+j] < m){
        G2_to_G1[j] = n; // epsilon -> j
      }
    }

    nged = this->GedFromMapping(g1, g2, G1_to_G2,n, G2_to_G1,m);

    if (ged > nged || ged == -1)
      ged = nged;
  }

#if XP_OUTPUT
  auto end = std::chrono::steady_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  std::cout << elapsed.count() << ":";
  std::cout << (mappings.size() == k) << ":";
#endif

  for (it=mappings.begin(); it!=mappings.end(); it++)
    delete[] *it;

  delete [] G1_to_G2;
  delete [] G2_to_G1;
  return ged;
}




 double RandomWalksGraphEditDistanceMulti::
 operator() (Graph<int,int> * g1,
             Graph<int,int> * g2)
{
  return operator() (g1,g2,nep);
}


RandomWalksGraphEditDistanceMulti::
~RandomWalksGraphEditDistanceMulti(){
  if (this->C != NULL)     delete [] this->C;
  if (this->Clsap != NULL) delete[] this->Clsap;
  this->C = NULL;
  this->Clsap = NULL;
}

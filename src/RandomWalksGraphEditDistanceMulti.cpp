/*
 * @file RandomWalksGraphEditDistanceMulti.cpp
 * @author Evariste <<evariste.daller@unicaen.fr>>
 * @version     0.0.1 - Wed Mar 29 2017
 */

 #include "RandomWalksGraphEditDistanceMulti.h"


void RandomWalksGraphEditDistanceMulti::
getOptimalMapping (Graph<int,int> * g1,
                   Graph<int,int> * g2,
                   int * G1_to_G2, int * G2_to_G1 )
{
  int n=g1->Size();
  int m=g2->Size();

  delete [] this->C;     //this->C = NULL;

#if XP_OUTPUT
  auto start = std::chrono::steady_clock::now();
#endif
  this->computeCostMatrix(g1, g2);

#if XP_OUTPUT
  auto end = std::chrono::steady_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  std::cout << elapsed.count() << ":";
#endif

  this->computeOptimalMapping(this, g1, g2, this->C, G1_to_G2, G2_to_G1);
}


double RandomWalksGraphEditDistanceMulti::
operator() (Graph<int,int> * g1,
            Graph<int,int> * g2)
{
  int n=g1->Size();
  int m=g2->Size();

  int* G1_to_G2 = new int[n];
  int* G2_to_G1 = new int[m];

  getOptimalMapping(g1, g2, G1_to_G2, G2_to_G1);

  delete [] G2_to_G1;
  delete [] G1_to_G2;

  return this->_ged;
}


RandomWalksGraphEditDistanceMulti::
~RandomWalksGraphEditDistanceMulti(){
  if (this->C != NULL)     delete [] this->C;
  this->C = NULL;
}

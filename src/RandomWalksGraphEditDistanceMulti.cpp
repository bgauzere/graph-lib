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
  clock_t t = clock();
#endif
  this->computeCostMatrix(g1, g2);

#if XP_OUTPUT
  t = clock() - t;
  _xp_out_ << ((float)t) / CLOCKS_PER_SEC << ":";
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

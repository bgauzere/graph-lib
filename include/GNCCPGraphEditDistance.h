/**
 * @file GNCCPGraphEditDistance.h
 * @author Benoit <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Wed Mar  8 2017
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 */

#ifndef __GNCCPGRAPHEDITDISTANCE_H__
#define __GNCCPGRAPHEDITDISTANCE_H__
#include <Eigen/Dense>
using namespace Eigen;

#include "GraphEditDistance.h"
#include "IPFPZetaGraphEditDistance.h"
#include "utils.h"

template<class NodeAttribute, class EdgeAttribute>
class GNCCPGraphEditDistance:
  public GraphEditDistance<NodeAttribute, EdgeAttribute>{
protected:
  
private:
  double _d = 0.1;
  double _zeta;
  IPFPZetaGraphEditDistance<NodeAttribute,EdgeAttribute> * sub_algo;
public:
  GNCCPGraphEditDistance(EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction):
    GraphEditDistance<NodeAttribute,EdgeAttribute>(costFunction){};

  virtual void getOptimalMapping(Graph<NodeAttribute,EdgeAttribute> * g1,
				 Graph<NodeAttribute,EdgeAttribute> * g2,
				 int * G1_to_G2, int * G2_to_G1);
  
  ~GNCCPGraphEditDistance(){}

};




template<class NodeAttribute, class EdgeAttribute>
void GNCCPGraphEditDistance<NodeAttribute, EdgeAttribute>::
getOptimalMapping(Graph<NodeAttribute,EdgeAttribute> * g1,
			       Graph<NodeAttribute,EdgeAttribute> * g2,
			       int * G1_to_G2, int * G2_to_G1){
  for(int i =0;i<g1->Size();i++)
    G1_to_G2[i] = i;
  for(int j=0;j<g2->Size();j++)
    G2_to_G1[j] = j;
  this->_zeta = 1;
  while(this->_zeta > -1){
    this->sub_algo = new IPFPZetaGraphEditDistance<NodeAttribute, EdgeAttribute>(this->cf, this->_zeta);
    this->sub_algo->IPFPalgorithm(g1,g2,G1_to_G2, G2_to_G1);
    this->_zeta -= this->_d;
    delete this->sub_algo; // pas top, faire evoluer zeta ?
  }
  
    
}





#endif // __GNCCPGRAPHEDITDISTANCE_H__

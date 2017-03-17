/**
 * @file BestPerfectMatching.cpp
 * @author Évariste <<21511575@etu.unicaen.fr>>
 */

#include "BestPerfectMatching.h"
#include <iostream>


BipartiteSCC::BipartiteSCC(unsigned int size_u, unsigned int size_v):
  u(size_u, false),
  v(size_v, false)
{}



void
BestPerfectMatching::strongConnect( const Eigen::MatrixXi& gm, int v ){

  vnum[v+offset] = num;
  vaccess[v+offset] = num;
  num += 1;

  tarjanStack.push(v);
  setstack.push(!(bool)(offset));
  instack[v+offset] = true;

  // For each w successor of v
  // 2 cases :
  // * v in X => 1 successor,  gm[v,w] ==  1
  // * v in Y => * successors, gm[w,v] == -1

  int w=0;
  if (offset == 0){ // v in X : find the 1
    while(w < gm.cols()){
      if (gm(v,w) == 1){
        if (vnum[w+offsize] == -1){ // if num(w) not defined
          // explore w
          offset = offsize; // w is in Y
          strongConnect(gm, w);
          offset = 0; // we are back in X
          vaccess[v] = std::min(vaccess[v], vaccess[w+offsize]);
        }
        else if(instack[w+offsize]){
          vaccess[v] = std::min(vaccess[v], vnum[w+offsize]);
        }
        //bypass the rest of search (there is max one successor)
        w = gm.cols();
      }
      w++;
    }
  }
  else{ // v in Y : find all the -1
    while(w < gm.rows()){
      if (gm(w,v) == -1){ // if w is a successor
        if (vnum[w] == -1){
          offset = 0; // w is in X
          strongConnect(gm, w);
          offset = offsize; // back in Y
          vaccess[v+offset] = std::min(vaccess[v+offset], vaccess[w]);
        }
        else if(instack[w]){
          vaccess[v+offset] = std::min(vaccess[v+offset], vnum[w]);
        }
      }
      w++;
    }
  }

  // If v is a root (v.access == v.num)
  if (vaccess[v+offset] == vnum[v+offset]){
    scc.emplace_back(gm.rows(), gm.cols());
    bool inX;
    do{
      if (!tarjanStack.empty()){
        w = tarjanStack.top();
        tarjanStack.pop();

        inX = setstack.top();
        setstack.pop();
        instack[w+offset] = false;

        if(inX)
          scc.back().u[w] = true;
        else
          scc.back().v[w] = true;
      }
    }while(w != v || !(bool)(offset) != inX);
  }
} // end strongConnect()


const std::vector<BipartiteSCC>&
BestPerfectMatching::findSCC( const Eigen::MatrixXi& gm ){

  // Initialization
  num = 0;
  access = 0;
  vnum.clear();
  vaccess.clear();
  instack.clear();

  scc.clear();

  vnum.resize(gm.rows()+gm.cols(), -1);
  instack.resize(gm.rows()+gm.cols(), false);
  vaccess.resize(gm.rows()+gm.cols());

  // Begin in X
  offset = 0;
  offsize = gm.rows();


  // For each vertice in X (rows of gm)
  for (int i=0; i<gm.rows(); i++){
    if (vnum[i] != -1)
      offset = 0;
      strongConnect(gm, i);
  }

  return scc;
}


void
BestPerfectMatching::rmUnnecessaryEdges( Eigen::MatrixXi& gm, const std::vector<BipartiteSCC>& scc_in ){
  for (int i=0; i<gm.rows(); i++){
    for (int j=0; j<gm.cols(); j++){
      bool inSomeSCC = false;
      for (int s=0; s<scc_in.size(); s++){
        inSomeSCC = inSomeSCC || (scc_in[s].u[i] && scc_in[s].v[j]);
      }
      gm(i,j) = gm(i,j) * inSomeSCC;
    }
  }
}


void
BestPerfectMatching::rmUnnecessaryEdges( Eigen::MatrixXi& gm){
  rmUnnecessaryEdges(gm, scc);
}

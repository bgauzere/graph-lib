/**
 * @file BestPerfectMatching.cpp
 * @author Évariste <<21511575@etu.unicaen.fr>>
 */

#include "BestPerfectMatching.h"



BipartiteSCC::BipartiteSCC(unsigned int size_u, unsigned int size_v):
  u(size_u, false),
  v(size_v, false)
{}



void
BestPerfectMatching::strongConnect( const Eigen::MatrixXi& gm, int v ){

  vnum[v+offset] = num;
  vaccess[v+offset] = num;
  num += 1;

  tarjanStack.push(v+offset);
  instack[v+offset] = true;

  // For each w successor of v
  // 2 cases :
  // * v in X => 1 successor,  gm[v,w] ==  1
  // * v in Y => * successors, gm[w,v] == -1

  int w=0;
  if (offset == 0){ // v in X : find the 1
    while(gm(v,w) != 1) w++;

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
    int offbakup = offset;
    do{
      offset = offsize-offset; //swap offset into 0 or offsize
      w = tarjanStack.top();
      tarjanStack.pop();
      if(offset)
        scc.back().v[w] = true;
      else
        scc.back().u[w] = true;
    }while(w != v || offset != offbakup);
  }
}


const std::vector<BipartiteSCC>&
BestPerfectMatching::findSCC( const Eigen::MatrixXi& gm ){

  // Initialization
  num = 0;
  access = 0;
  vnum.clear();
  vaccess.clear();
  instack.clear();

  scc.clear();
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
      strongConnect(gm, i);
  }

  return scc;
}


void
BestPerfectMatching::rmUnnecessaryEdges( Eigen::MatrixXi& gm, const std::vector<BipartiteSCC>& scc ){
  //TODO
}


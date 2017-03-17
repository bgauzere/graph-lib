#include <iostream>
#include <Eigen/Dense>
#include "BestPerfectMatching.h"

using namespace std;
using namespace Eigen;

int main(){

  /************* Exemple 1 *******************/

  MatrixXi gm(4,4);
  gm << 0,  1,  0, -1,
       -1,  0,  0,  0,
        0, -1,  1, -1,
        0,  0, -1,  1;

  cout << gm << endl;

  BestPerfectMatching bpm;
  const vector<BipartiteSCC>& res = bpm.findSCC(gm);

  cout << "found " << res.size() << " scc" << endl;

  for (int i=0; i<res.size(); i++){
    cout << "u_" << i << " : ";
    for (int j=0; j<res[i].u.size(); j++){
      cout << res[i].u[j] << " ";
    }
    cout << endl;

    cout << "v_" << i << " : ";
    for (int j=0; j<res[i].v.size(); j++){
      cout << res[i].v[j] << " ";
    }
    cout << endl;
  }

  cout << endl;

  /*** Test rmUnnecessaryEdges ***/

  //bpm.rmUnnecessaryEdges(gm, res);
  bpm.rmUnnecessaryEdges(gm);

  cout << endl;
  cout << "After rmUnnecessaryEdges : " << endl;
  cout << gm << endl;

  cout << endl << "-----" << endl;

  /************* Exemple 2 *******************/

  MatrixXi am(3,3);
  am << 1,  0,  0,
       -1, -1,  1,
        0,  1, -1;

  cout << am << endl;
  const std::vector<BipartiteSCC>& res2 = bpm.findSCC(am);

  cout << "found " << res2.size() << " scc" << endl;

  for (int i=0; i<res2.size(); i++){
    cout << "u_" << i << " : ";
    for (int j=0; j<res2[i].u.size(); j++){
      cout << res2[i].u[j] << " ";
    }
    cout << endl;

    cout << "v_" << i << " : ";
    for (int j=0; j<res2[i].v.size(); j++){
      cout << res2[i].v[j] << " ";
    }
    cout << endl;
  }

  //bpm.rmUnnecessaryEdges(gm, res);
  bpm.rmUnnecessaryEdges(am);

  cout << endl;
  cout << "After rmUnnecessaryEdges : " << endl;
  cout << am << endl;

  cout << endl << "-----" << endl;

  /************* Exemple 3 *******************/

  MatrixXi bm(5,4);
  bm << 1, -1,  0,  0,
        0, -1,  1, -1,
       -1,  1,  0,  0,
        0,  0, -1,  1,
        0,  0,  0, -1;

  cout << bm << endl;
  const std::vector<BipartiteSCC>& res3 = bpm.findSCC(bm);

  cout << "found " << res3.size() << " scc" << endl;

  for (int i=0; i<res3.size(); i++){
    cout << "u_" << i << " : ";
    for (int j=0; j<res3[i].u.size(); j++){
      cout << res3[i].u[j] << " ";
    }
    cout << endl;

    cout << "v_" << i << " : ";
    for (int j=0; j<res3[i].v.size(); j++){
      cout << res3[i].v[j] << " ";
    }
    cout << endl;
  }

  //bpm.rmUnnecessaryEdges(gm, res);
  bpm.rmUnnecessaryEdges(bm);

  cout << endl;
  cout << "After rmUnnecessaryEdges : " << endl;
  cout << bm << endl;
  return 0;
}

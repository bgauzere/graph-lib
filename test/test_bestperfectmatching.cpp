#include <iostream>
#include <Eigen/Dense>
#include "../include/BestPerfectMatching.h"

using namespace std;
using namespace Eigen;

int main(){

  MatrixXi gm(4,4);
  gm << 0, 1, 0, -1,
       -1, 0, 0, 0,
       0, -1, 1, -1,
       0, 0, -1, 1;

  cout << gm << endl;

  const vector<BipartiteSCC>& res = BestPerfectMatching::findSCC(gm);

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

  return 0;
}

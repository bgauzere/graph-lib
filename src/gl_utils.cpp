/*
 * @file utils.cpp
 * @author Benoit <<benoit.gauzere@insa-rouen.fr>>
 * @version     0.0.1 - Thu Feb  2 2017
 *
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *
 * Description of the program objectives.
 * All necessary references.
 *
 */
#include <cstring>
#include <iostream>

#include "gl_utils.h"

std::vector<char*> split (const char* chaine, const char* sep){
  std::vector<char*> v;

  char* s = strtok ((char*)chaine, (char*)sep);

  while (s != NULL)
    {
      v.push_back (s);
      s = strtok (NULL, (char*)sep);
    }

  return v;
}

//int sub2ind(int i, int j, int n){return i + j*n;};

void print_array(double* arr, int n, int m){
  for (int i=0; i<n; i++){
    for (int j=0; j<m; j++)
      std::cout << arr[i + j*n] << "  ";
    std::cout << std::endl;
  }
}

void mappings_lsap2lsape( unsigned int* lsapMapping,
			  unsigned int* rho,
			  unsigned int* varrho,
			  int n, int m ){
  for (int j=0; j<m; j++){ // connect all to epsilon by default
    varrho[j] = n;
  }

  for (int i=0; i<n; i++){
    if ((int)(lsapMapping[i]) >= m)
      rho[i] = m; // i -> epsilon
    else{
      rho[i] = lsapMapping[i];
      varrho[lsapMapping[i]] = i;
    }
  }
}



void mappings_lsap2lsape( std::list<unsigned int*> lsapMaps,
			  std::list<unsigned int*>& G1toG2,
			  std::list<unsigned int*>& G2toG1,
			  int n, int m ){
  G1toG2.clear();
  G2toG1.clear();
  for ( std::list<unsigned int*>::iterator it=lsapMaps.begin();
        it != lsapMaps.end();  it++ )
  {
    unsigned int* lsapMapping = *it;
    unsigned int* rho = new unsigned int[n];
    unsigned int* varrho = new unsigned int[m];

    mappings_lsap2lsape(lsapMapping, rho, varrho, n, m);

    G1toG2.push_back(rho);
    G2toG1.push_back(varrho);
  }
}


//double abs(double x) { return (x >= 0) ? x : -x; }

/**
 * @file gl_utils.h
 * @author Benoit <<benoit.gauzere@insa-rouen.fr>> 
 * @version     0.0.1 - Thu Feb  2 2017
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 */

#ifndef __GL_UTILS_H__
#define __GL_UTILS_H__

#include <vector>
#include <list>

inline bool is_mapping_valid(unsigned int mapping, unsigned int nb_nodes){
  return (mapping < nb_nodes);
}

inline int  sub2ind(int i, int j, int n){ return  (i + (j) * (n));}

std::vector<char*> split (const char* chaine, const char* sep);

//int sub2ind(int i, int j, int n);

template<typename T>
T mean(T * tab, int size){
  T sum = 0;
  for (int i =0;i<size;i++)
    sum += tab[i];
  return sum/size;
}


void mappings_lsap2lsape( unsigned int* lsapMapping,
			  unsigned int* rho,
			  unsigned int* varrho,
			  int n, int m );

void mappings_lsap2lsape( std::list<unsigned int*> lsapMaps,
			  std::list<unsigned int*>& G1toG2,
			  std::list<unsigned int*>& G2toG1,
			  int n, int m );

#endif // __UTILS_H__

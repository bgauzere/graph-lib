/*
 * @file utils.cpp
 * @author Benoit <<benoit.gauzere@greyc.ensicaen.fr>> 
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

#include "utils.h"

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

int sub2ind(int i, int j, int n){return i + j*n;};


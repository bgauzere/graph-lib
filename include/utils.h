/**
 * @file utils.h
 * @author Benoit <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Thu Feb  2 2017
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 */

#ifndef __UTILS_H__
#define __UTILS_H__

#include <vector>

std::vector<char*> split (const char* chaine, const char* sep);
//TODO:a inliner
int sub2ind(int i, int j, int n);

#endif // __UTILS_H__

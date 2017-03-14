/**
 * @file SymbolicGraph.h
 * @author Benoit <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Sat Feb  4 2017
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 */

#ifndef __SYMBOLICGRAPH_H__
#define __SYMBOLICGRAPH_H__
#include "graph.h"
class SymbolicGraph: public Graph<int,int>
{
private:
  static int readChemicalEdgeLabel(TiXmlElement *elem);
  static int readChemicalNodeLabel(TiXmlElement *elem);


public:
  SymbolicGraph(const char * filename);
  SymbolicGraph(int * am, int nb_nodes, bool directed);
  int * getLabeledAdjacencyMatrix();
};

#endif // __SYMBOLICGRAPH_H__


/**
 * @file Dataset.h
 * @author Benoit <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Mon Feb  6 2017
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 */

#ifndef __DATASET_H__
#define __DATASET_H__
#include <vector>
#include <string>
#include <libgen.h>
#include "graph.h"
#include "SymbolicGraph.h"
#include "GraphEditDistance.h"

template <class NodeAttribute, class EdgeAttribute, class PropertyType>
class Dataset
{
private:
  std::vector<PropertyType> properties;
  std::vector<Graph<NodeAttribute, EdgeAttribute> * > graphs;
  // void loadCXL(const char * filename);
public:
  //  Dataset(const char * filename){};
  Graph<NodeAttribute, EdgeAttribute> * getGraph(unsigned int id){return graphs[id];};
  Graph<NodeAttribute, EdgeAttribute> * operator[](unsigned int id){return graphs[id];};

  PropertyType getProperty(unsigned int id){return properties[id];};
  PropertyType operator()(unsigned int id){return properties[id];};

  int add(Graph<NodeAttribute, EdgeAttribute> * g, PropertyType y){graphs.push_back(g);properties.push_back(y);return graphs.size();}

  int size(){return graphs.size();};

  double * computeGraphEditDistance(GraphEditDistance<NodeAttribute,EdgeAttribute> * ed);
  void shuffleize();
};


template<class NodeAttribute,class EdgeAttribute, class PropertyType>
void Dataset<NodeAttribute,EdgeAttribute,PropertyType>::shuffleize(){
  int N = size();
  for (int i=0;i<N;i++)
    (*this)[i]->shuffleize();
}

template<class NodeAttribute,class EdgeAttribute, class PropertyType>
double * Dataset<NodeAttribute,EdgeAttribute,PropertyType>::computeGraphEditDistance(GraphEditDistance<NodeAttribute,EdgeAttribute> * ed){
  int N = size();
  double * distances = new double[size()*size()];
  for (int i=0;i<N;i++)
    for (int j=0;j<N;j++)
      distances[sub2ind(i,j,N)] = (*ed)(graphs[i],graphs[j]);

  return distances;
}




template<class PropertyType>
class ChemicalDataset: public Dataset<int,int,PropertyType>
{
private:
  void loadDS(const char * filename);
public:
  ChemicalDataset(const char * filename);
};
  
template<class PropertyType>
void ChemicalDataset<PropertyType>::loadDS(const char* filename){
  
  std::ifstream f_tmp (filename);
  char * unconst_filename = new char[strlen(filename)+1];
  unconst_filename = strcpy(unconst_filename, filename);
  char * path = dirname(unconst_filename);
  if (f_tmp.is_open()){
    std::string s;
    while (getline(f_tmp, s))
      if (s[0] != '#'){
	std::string path_ctfile(path);
	path_ctfile += std::string("/");
	std::istringstream liness(s);
	std::string ctfile;
	liness >> ctfile;
	PropertyType y;
	liness >> y;
	std::string full_ctfile = path_ctfile;
	full_ctfile += ctfile;
	SymbolicGraph * g = new SymbolicGraph(full_ctfile.c_str());
	this->add(g,y);
      }
  }
  f_tmp.close();
  
}
template<class PropertyType>
ChemicalDataset<PropertyType>::ChemicalDataset(const char * filename){
  const char * ext = strrchr(filename,'.'); 
  if (strcmp(ext,".ds") == 0){
    loadDS(filename);
  }
}


#endif // __DATASET_H__

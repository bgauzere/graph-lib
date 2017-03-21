/*
 * @file test_graph.cpp
 * @author Benoit <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Wed Feb  1 2017
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 *
 */

#include <iostream>
#include "graph.h"
#include "tinyxml.h"
#include "Dataset.h"
#include "GraphEditDistance.h"
#include "SymbolicGraph.h"
#include "ConstantGraphEditDistance.h"
#include "BipartiteGraphEditDistance.h"
#include "RandomWalksGraphEditDistance.h"
#include "IPFPGraphEditDistance.h"
#include "GNCCPGraphEditDistance.h"

using namespace std;

typedef struct{
  int x;
  int y;
  string type;
} GrecNodeLabel;

typedef struct{
  int frequency;
  string type;
} GrecEdgeLabel;


GrecEdgeLabel * readGrecEdgeLabel(TiXmlElement *elem){
  GrecEdgeLabel * result = new GrecEdgeLabel;
  TiXmlElement* child = elem->FirstChildElement();
  while ( child ){
    string childName =  child->Attribute("name");
    TiXmlElement * child2 = child->FirstChildElement();
    if ( child2 ){
      if(childName.compare("frequency")==0){
	result->frequency = std::stof(child2->GetText());
      }
      if(childName.compare("type0")==0){
	result->type = child2->GetText();
      }
    }
    child = child->NextSiblingElement(); // iteration
  }
  return result;
}

GrecNodeLabel * readGrecNodeLabel(TiXmlElement *elem){
  GrecNodeLabel * result = new GrecNodeLabel;
  TiXmlElement* child = elem->FirstChildElement();
  while ( child ){
    string childName =  child->Attribute("name");
    TiXmlElement * child2 = child->FirstChildElement();
    if ( child2 ){
      if(childName.compare("x")==0){
	result->x = stod(child2->GetText());
      }else if (childName.compare("y")==0){
	result->y = stod(child2->GetText());
      } else if (childName.compare("type")==0){
	result->type = child2->GetText();;
      }
    }
    child = child->NextSiblingElement(); // iteration au prochain attribut
    
  }
  return result;
}


void usage (char * s)
{
  cerr << "Usage : "<< s << endl;
  cerr << "options:" << endl;
}



int main (int argc, char** argv)
{

  
  ConstantEditDistanceCost * cf = new ConstantEditDistanceCost(1,3,3,1,3,3);
  GNCCPGraphEditDistance<int,int> *ed_gnccp = new GNCCPGraphEditDistance<int,int>(cf);

  BipartiteGraphEditDistance<int,int> *ed_bunke = new BipartiteGraphEditDistance<int,int>(cf);
  RandomWalksGraphEditDistance *ed_rw = new RandomWalksGraphEditDistance(cf,3);

  IPFPGraphEditDistance<int,int> *ed_ipfpe_rw = new IPFPGraphEditDistance<int,int>(cf,ed_rw);
  IPFPGraphEditDistance<int,int> *ed_ipfpe_bunke = new IPFPGraphEditDistance<int,int>(cf,ed_bunke);

  // ChemicalDataset<double> AcyclicDataset("/home/bgauzere/work/Datasets/Acyclic/dataset_bps.ds");
  // AcyclicDataset.shuffleize();
  // double * distancesAcyclic_IPFPE_RW =  AcyclicDataset.computeGraphEditDistance(ed_ipfpe_rw);
  // double * distancesAcyclic_IPFPE_BUNKE =  AcyclicDataset.computeGraphEditDistance(ed_ipfpe_bunke);
  // double * distancesAcyclic_LSAPE_RW =  AcyclicDataset.computeGraphEditDistance(ed);
  // double * distancesAcyclic_LSAPE_BUNKE =  AcyclicDataset.computeGraphEditDistance(ed_rw);
  // cout << AcyclicDataset.size()<< endl;
  // ofstream f("gram.mat");
  // double sum_bunke =0.0,
  //   sum_rw=0.0,sum_lsape_rw=0.0,sum_lsape_bunke=0.0;
  // for(int i=0;i<AcyclicDataset.size();i++){
  //   for(int j=0;j<AcyclicDataset.size();j++){
  //     sum_bunke += distancesAcyclic_IPFPE_BUNKE[sub2ind(i,j,AcyclicDataset.size())];
  //     sum_rw += distancesAcyclic_IPFPE_RW[sub2ind(i,j,AcyclicDataset.size())];
  //     sum_lsape_bunke += distancesAcyclic_LSAPE_BUNKE[sub2ind(i,j,AcyclicDataset.size())];
  //     sum_lsape_rw += distancesAcyclic_LSAPE_RW[sub2ind(i,j,AcyclicDataset.size())];
  //     f<< distancesAcyclic_IPFPE_RW[sub2ind(i,j,AcyclicDataset.size())];
  //     if(j != AcyclicDataset.size() -1)
  //      	f << ",";
  //   }
  //   f << endl;
  // }
  // f.close();
  
  
  // cout << "Average distance pour LSAPE Bunke sur Acyclic : "<< sum_lsape_bunke/(AcyclicDataset.size()*AcyclicDataset.size()) <<endl;
  // cout << "Average distance pour LSAPE RW sur Acyclic : "<< sum_lsape_rw/(AcyclicDataset.size()*AcyclicDataset.size()) <<endl;
  // cout << "Average distance pour IPFPE Bunke sur Acyclic : "<< sum_bunke/(AcyclicDataset.size()*AcyclicDataset.size()) <<endl;
  // cout << "Average distance pour IPFPE RW sur Acyclic : "<< sum_rw/(AcyclicDataset.size()*AcyclicDataset.size()) <<endl;
  

  // Graph<GrecNodeLabel*, GrecEdgeLabel*> *g2 = new Graph<GrecNodeLabel*,GrecEdgeLabel *>("./grec.gxl", readGrecNodeLabel,readGrecEdgeLabel);
  // cout << "Must be 12 : " << g2->Size()  << endl;
  // cout << (*g2)[0]->attr->x << endl;
  // cout << (*g2)[11]->attr->type << endl;

  
  // SymbolicGraph *g3 = new SymbolicGraph("test.ct");
  // cout << "Must be 9 : " << g3->Size()  << endl;
  
  // SymbolicGraph *g4 = new SymbolicGraph("molecule.gxl");
  // cout << "Must be ? : " << g4->Size()  << endl;

  
  

  // ConstantDataset<double> mydataset("./dataset.ds");
  // cout << mydataset.size() << endl;
  // cout << mydataset[1]->Size() << endl;
  // cout << mydataset(1) << endl;


  
  int g5_am[9] = {6,0,1,
  		  0,6,1,
  		  1,1,8};
  int g5_size = 3;
  SymbolicGraph * g5 = new SymbolicGraph((int*) g5_am, g5_size,false);
  
  int g6_am[16] = {6,0,1,0,
  		   0,6,0,1,
  		   1,0,8,1,
  		   0,1,1,8};
  int g6_size = 4;
  SymbolicGraph * g6 = new SymbolicGraph((int *)g6_am, g6_size,false);

  // cout << "Size de g6" << g6->Size() << endl;
  // cout << "Nb Edges de g6" << g6->getNbEdges() << endl;

  // cout << g5->Size() << endl;
  // cout << g5->getNbEdges() << endl;

  // int g5tog6 [3] = {0,1,2};
  // int g6tog5 [4] = {0,1,2,3};
  
  
  // cout << "Distance between g5 and g6 " << ed->GedFromMapping(g5,g6,(int *)g5tog6,g5->Size(),(int*)g6tog5, g6_size) << endl;
  // cout << "Distance between g6 and g5 " << ed->GedFromMapping(g6,g5,(int *)g6tog5,g6->Size(),(int*)g5tog6,g5_size) << endl;
 
  // int g5tog5 [3] = {0,1,2};
  // cout << "Distance between g5 and itself :" << ed->GedFromMapping(g5,g5,(int *)g5tog5, g5_size, (int *)g5tog5, g5_size ) << endl;

  // int g6tog6 [8] = {0,1,2,3};
  // cout << "Distance between g6 and itself :" << ed->GedFromMapping(g6,g6,(int *)g6tog6, g6_size,(int *)g6tog6, g6_size) << endl;
  
  // int * AM = g6->getLabeledAdjacencyMatrix();
  // for (int i =0;i<4;i++){
  //   for (int j =0;j<4;j++)
  //     cout << AM[sub2ind(i,j,4)]<< ",";
  //   cout << endl;
  // }
  // //g6->shuffleize();
  // cout << endl;
  // AM = g6->getLabeledAdjacencyMatrix();
  // for (int i =0;i<4;i++){
  //   for (int j =0;j<4;j++)
  //     cout << AM[sub2ind(i,j,4)]<< ",";
  //   cout << endl;
  // }
  //g5->shuffleize();
  //cout << "Distance between g6 and g5 avec RW : " << (*ed_rw)(g5,g6)  << endl;


  cout << "Distance between g5 and g6 avec LSAPE Bunke : " << (*ed_bunke)(g5,g6)  << endl;
  cout << "Distance between g5 and g6 avec LSAPE RW    : " << (*ed_rw)(g5,g6)  << endl;
  cout << "Distance between g5 and g6 avec IPFPE Bunke : " << (*ed_ipfpe_bunke)(g5,g6)  << endl;
  cout << "Distance between g5 and g6 avec IPFPE RW    : " << (*ed_ipfpe_rw)(g5,g6)  << endl;
  cout << "Distance between g5 and g6 avec GNCCP       : " << (*ed_gnccp)(g5,g6)  << endl;

  // cout << "Distance between g5 and g6 avec IPFPE : " << (*ed_ipfpe)(g5,g6)  << endl;
  
  
  // cout << "Average distance (min) pour LSAP Bunke sur Acyclic : "<< sum_bis/(AcyclicDataset.size()*AcyclicDataset.size()) <<endl;

  // distancesAcyclic =  AcyclicDataset.computeGraphEditDistance(ed_rw);
  // sum =0.0; sum_bis=0.0;
  // for(int i=0;i<AcyclicDataset.size();i++){
  //   for(int j=0;j<AcyclicDataset.size();j++){
  //     sum += distancesAcyclic[sub2ind(i,j,AcyclicDataset.size())];
  //     sum_bis += min(distancesAcyclic[sub2ind(i,j,AcyclicDataset.size())],distancesAcyclic[sub2ind(j,i,AcyclicDataset.size())]);
  //     // f<< distancesAcyclic[sub2ind(i,j,AcyclicDataset.size())];
  //     // if(j != AcyclicDataset.size() -1)
  //     // 	f << ",";
  //   }
  //   // f << endl;
  // }
  
  // cout << "Average distance pour LSAP RW sur Acyclic : "<< sum/(AcyclicDataset.size()*AcyclicDataset.size()) <<endl;
  
  // cout << "Average distance (min) pour LSAP RW sur Acyclic : "<< sum_bis/(AcyclicDataset.size()*AcyclicDataset.size()) <<endl;

  
  // ChemicalDataset<double> PAHDataset("/home/bgauzere/work/Datasets/PAH/dataset.ds");
  // double * distancesPAH =  PAHDataset.computeGraphEditDistance(ed);
  // cout << PAHDataset.size()<< endl;
  // sum =0.0;
  // sum_bis=0.0;
  // for(int i=0;i<PAHDataset.size();i++){
  //   for(int j=0;j<PAHDataset.size();j++){
  //     sum += distancesPAH[sub2ind(i,j,PAHDataset.size())];
  //     sum_bis += min(distancesPAH[sub2ind(i,j,PAHDataset.size())],distancesPAH[sub2ind(j,i,PAHDataset.size())]);
  //   }
    
  // }
  // cout << "Average distance pour LSAP Bunke sur PAH : "<< sum/(PAHDataset.size()*PAHDataset.size()) <<endl;
  // cout << "Average distance (min) pour LSAP Bunke sur PAH : "<< sum_bis/(PAHDataset.size()*PAHDataset.size()) <<endl;
  delete ed_bunke;
  delete ed_rw;
  delete ed_ipfpe_bunke;
  delete ed_ipfpe_rw;
  delete ed_gnccp;
  delete cf;
  delete g6;
  delete g5;
  return 0;
}

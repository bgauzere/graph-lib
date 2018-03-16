#include <unistd.h>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>


#include "graph.h"
#include "tinyxml.h"
#include "Dataset.h"
#include "GraphEditDistance.h"
#include "SymbolicGraph.h"
#include "ConstantGraphEditDistance.h"
#include "BipartiteLowerBoundGraphEditDistance.h"
#include "BipartiteGraphEditDistance.h"
#include "gl_utils.h"

#include <lsape.h>

using namespace std;


void usage (char * s)
{
  cerr << "Usage : "<< s << " dataset " << " options"<<endl;
  cerr << "options:" << endl;
  cerr << "\t -m solver " << endl;
  cerr << "\t \t Specify the solver to resolve LSAP" << endl;
  cerr << "\t -o output_file " << endl;
  cerr << "\t \t Specify a filename for outputing lower bounds matrix" << endl;
  cerr << "\t -c cns,cni,cnd,ces,cei,ced " << endl;
  cerr << "\t \t Specify edit operation costs" << endl;
}

struct Options{
  string dataset_file = "";
  string method = "";
  string output_file = "";
  double cns = 1;
  double cni = 3;
  double cnd = 3;
  double ces = 1;
  double cei = 3;
  double ced = 3;
  bool shuffle = false;
};

void saveMatrix(double * matrix,int N, const char * filename){
  std::filebuf fb;
  fb.open (filename,std::ios::out);
  std::ostream os(&fb);
  for (int i=0; i<N; i++){
    for (int j=0; j<N; j++){
      os << matrix[sub2ind(i,j,N)];
      if(j < N-1)
	os << ",";
    }
    os << endl;
  }
  fb.close();
}

struct Options * parseOptions(int argc, char** argv){
  struct Options * options = new struct Options();
  options->dataset_file = string(argv[1]);
  int opt;
  stringstream sstream;
  while ((opt = getopt(argc, argv, "m:o:c:sp:z")) != -1) {
    switch (opt) {
    case 'm':
      options->method = string(optarg);
      break;
    case 'o':
      options->output_file = string(optarg);
      break;
    case 'c':
      char tmp;
      sstream << optarg;
      sstream >> options->cns;      sstream >> tmp;
      sstream >> options->cni;      sstream >> tmp;
      sstream >> options->cnd;      sstream >> tmp;
      sstream >> options->ces;      sstream >> tmp;
      sstream >> options->cei;      sstream >> tmp;
      sstream >> options->ced;
      break;
    case 's':
      options->shuffle=true;
      break;
    default: /* '?' */
      cerr << "Options parsing failed."  << endl;
      usage(argv[0]);
      exit(EXIT_FAILURE);
    }
  }
  return options;
}

template <class NodeAttribute, class EdgeAttribute, class PropertyType>
void computeAndSaveUpperBounds(Dataset< NodeAttribute, EdgeAttribute, PropertyType> * dataset,
			       BipartiteGraphEditDistance<NodeAttribute, EdgeAttribute> * ub,
			       string prefix_filename){

  int N = dataset->size();
  double* upper_bounds = new double[N*N];
  double* times = new double[N*N];
  for (int i=0; i<N; i++){
    for (int j=i; j<N; j++){
      cout << "Calcul de lb(" << i << "," << j << ")"<< endl;
      Graph<NodeAttribute,EdgeAttribute> * g1 = (*dataset)[i];
      Graph<NodeAttribute,EdgeAttribute> * g2 = (*dataset)[j];
      double compTime;
      clock_t t = clock();
      if(g1->Size() > g2->Size())
	upper_bounds[sub2ind(j,i,N)]=upper_bounds[sub2ind(i,j,N)] = ub->getUpperBound(g2,g1,compTime);
      else
	upper_bounds[sub2ind(j,i,N)]=upper_bounds[sub2ind(i,j,N)] = ub->getUpperBound(g1,g2,compTime);

      times[sub2ind(j,i,N)]=times[sub2ind(i,j,N)] = compTime;
      
    }
  }
  saveMatrix(upper_bounds,N, (prefix_filename + string("ub.mat")).c_str());
  saveMatrix(times,N, (prefix_filename + string("ub_times.mat")).c_str());
}

template <class NodeAttribute, class EdgeAttribute, class PropertyType>
void computeAndSaveLowerBounds(Dataset< NodeAttribute, EdgeAttribute, PropertyType> * dataset,
			       BipartiteLowerBoundEditDistance<NodeAttribute, EdgeAttribute> * lb,
			       string prefix_filename){

  int N = dataset->size();
  double* lower_bounds = new double[N*N];
  double* times = new double[N*N];
  
  // Graph<NodeAttribute,EdgeAttribute> * g1 = (*dataset)[0];
  // Graph<NodeAttribute,EdgeAttribute> * g2 = (*dataset)[0];
  // double compTime;
  // lb->getLowerBound(g1,g2,compTime);
  // return;
  
  for (int i=0; i<N; i++){
    for (int j=i; j<N; j++){
      cout << "Calcul de lb(" << i << "," << j << ")"<< endl;
      Graph<NodeAttribute,EdgeAttribute> * g1 = (*dataset)[i];
      Graph<NodeAttribute,EdgeAttribute> * g2 = (*dataset)[j];
      double compTime;
      if(g1->Size() > g2->Size())
  	lower_bounds[sub2ind(j,i,N)]=lower_bounds[sub2ind(i,j,N)] = lb->getLowerBound(g2,g1,compTime);
      else
  	lower_bounds[sub2ind(j,i,N)]=lower_bounds[sub2ind(i,j,N)] = lb->getLowerBound(g1,g2,compTime);

      times[sub2ind(j,i,N)]=times[sub2ind(i,j,N)] = compTime;
      // cout << (int)lower_bounds[sub2ind(i,j,N)];
      // cout << endl;
    }
  }
  saveMatrix(lower_bounds,N, (prefix_filename + string("lb.mat")).c_str());
  saveMatrix(times,N, (prefix_filename + string("lb_times.mat")).c_str());

}



int main (int argc, char** argv)
{

  struct Options * options =   parseOptions(argc,argv);

  ConstantEditDistanceCost* cf = new ConstantEditDistanceCost(options->cns,options->cni, options->cnd,
							       options->ces,options->cei, options->ced);

  
  ChemicalDataset<double> * dataset = new ChemicalDataset<double>(options->dataset_file.c_str());

 
  if(options->shuffle)
    dataset->shuffleize();

  computeAndSaveUpperBounds(dataset,new BipartiteGraphEditDistance<int,int>(cf,lsape::ECBP),
  			    options->output_file+string("_ECBP_"));
  computeAndSaveUpperBounds(dataset,new BipartiteGraphEditDistance<int,int>(cf,lsape::FLWC),
  			    options->output_file+string("_FLWC_"));
  computeAndSaveUpperBounds(dataset,new BipartiteGraphEditDistance<int,int>(cf,lsape::EBP),
  			    options->output_file+string("_EBP_"));
  computeAndSaveUpperBounds(dataset,new BipartiteGraphEditDistance<int,int>(cf,lsape::FLCC),
  			    options->output_file+string("_FLCC_"));
  computeAndSaveUpperBounds(dataset,new BipartiteGraphEditDistance<int,int>(cf,lsape::FBP),
  			    options->output_file+string("_FBP_"));
  computeAndSaveUpperBounds(dataset,new BipartiteGraphEditDistance<int,int>(cf,lsape::FBP0),
  			    options->output_file+string("_FBP0_"));
  computeAndSaveUpperBounds(dataset,new BipartiteGraphEditDistance<int,int>(cf,lsape::SFBP),
  			    options->output_file+string("_SFBP_"));
  
  computeAndSaveLowerBounds(dataset,new BipartiteLowerBoundEditDistance<int,int>(cf,lsape::ECBP),
  			    options->output_file+string("_ECBP_"));
  computeAndSaveLowerBounds(dataset,new BipartiteLowerBoundEditDistance<int,int>(cf,lsape::FLWC),
  			    options->output_file+string("_FLWC_"));
  computeAndSaveLowerBounds(dataset,new BipartiteLowerBoundEditDistance<int,int>(cf,lsape::EBP),
  			    options->output_file+string("_EBP_"));
  computeAndSaveLowerBounds(dataset,new BipartiteLowerBoundEditDistance<int,int>(cf,lsape::FLCC),
  			    options->output_file+string("_FLCC_"));
  computeAndSaveLowerBounds(dataset,new BipartiteLowerBoundEditDistance<int,int>(cf,lsape::FBP),
  			    options->output_file+string("_FBP_"));
  computeAndSaveLowerBounds(dataset,new BipartiteLowerBoundEditDistance<int,int>(cf,lsape::FBP0),
  			    options->output_file+string("_FBP0_"));
  computeAndSaveLowerBounds(dataset,new BipartiteLowerBoundEditDistance<int,int>(cf,lsape::SFBP),
  			    options->output_file+string("_SFBP_"));
 
  delete dataset;
  delete cf;
  delete options;

  return 0;
}

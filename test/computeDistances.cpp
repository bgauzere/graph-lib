/*
 * @file test_graph.cpp
 * @author Évariste <<evariste.daller@unicaen.fr>>
 *
 *
 * Calcule toutes les distances entre toutes les molécules d'un datatset
 *
 */


#include <unistd.h>
#include <iostream>
#include <string>
#include <sstream>


#include "graph.h"
#include "tinyxml.h"
#include "Dataset.h"
#include "GraphEditDistance.h"
#include "SymbolicGraph.h"
#include "ConstantGraphEditDistance.h"
#include "CMUCostFunction.h"
#include "BipartiteGraphEditDistance.h"
#include "BipartiteGraphEditDistanceMulti.h"
#include "GreedyGraphEditDistance.h"
#include "RandomWalksGraphEditDistance.h"
#include "RandomWalksGraphEditDistanceMulti.h"
#include "IPFPGraphEditDistance.h"
#include "RandomMappings.h"
#include "MultistartRefinementGraphEditDistance.h"
#include "GNCCPGraphEditDistance.h"
#include "gl_utils.h"
using namespace std;



void usage (char * s)
{
  cout << "Usage : "<< s << " dataset   [options]  -m  method " << endl;
  cout << "options:" << endl;
  cout << "\t -m method \t \t Specify the algorithm used to compute edit distance" << endl;
  cout << "\t -c cns,cni,cnd,ces,cei,ced\t Specify edit operation costs" << endl;
  cout << "\t -p n_edit_paths \t \t Specify the number of edit paths to compute GED (lsape_multi)" << endl;
  cout << "\t -s        \t \t Shuffle the graphs" << endl;
  cout << endl;
  cout << " method can be one of the following :" << endl;
  cout << "\t lsape_bunke" << endl;
  cout << "\t lsape_multi_bunke" << endl;
  cout << "\t lsape_rw" << endl;
  cout << "\t lsape_multi_rw" << endl;
  cout << "\t lsape_multi_greedy" << endl << endl;
  cout << "\t ipfpe_flat" << endl;
  cout << "\t ipfpe_bunke" << endl;
  cout << "\t ipfpe_multi_bunke " << endl;
  cout << "\t ipfpe_rw" << endl;
  cout << "\t ipfpe_multi_rw" << endl;
  cout << "\t ipfpe_multi_greedy" << endl;
  cout << "\t ipfpe_multi_random" << endl;
  cout << "\t gnccp" << endl;
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
  bool cmu = false;
  int k = 3;
  int nep = 100; // number of edit paths for allsolution
};

struct Options * parseOptions(int argc, char** argv){
  if (argc < 4){
    usage(argv[0]);
    exit(0);
  }

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
      sstream << optarg;
      sstream >> options->cns;
      sstream >> options->cni;
      sstream >> options->cnd;
      sstream >> options->ces;
      sstream >> options->cei;
      sstream >> options->ced;
      break;
    case 's':
      options->shuffle=true;
      break;
    case 'p':
      sstream << optarg;
      sstream >> options->nep;
      break;
    case 'z':
      options->cmu = true;
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
double * computeGraphEditDistance(Dataset< NodeAttribute, EdgeAttribute, PropertyType> * dataset,
				  GraphEditDistance<NodeAttribute, EdgeAttribute> * ed,
				  bool shuffle, int nep=0){
  if(shuffle)
    dataset->shuffleize();

  //double * distances =  dataset->computeGraphEditDistance(ed, true);
  int N = dataset->size();
  double* distances = new double[N*N];
  //struct timeval  tv1, tv2;
  for (int i=0; i<N; i++){
    for (int j=0; j<N; j++){
      #ifdef PRINT_TIMES
        gettimeofday(&tv1, NULL);
      #endif

      distances[sub2ind(i,j,N)] = (*ed)((*dataset)[i], (*dataset)[j]);;

      #ifdef PRINT_TIMES
        gettimeofday(&tv2, NULL);
        cout << ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec)) << ", " ;
      #endif

      cout << (int)distances[sub2ind(i,j,N)];
      cout << endl;
    }
  }
  return distances;
}



int main (int argc, char** argv)
{

  struct Options * options =   parseOptions(argc,argv);


  options->k = 3;
  ConstantEditDistanceCost* cf = new ConstantEditDistanceCost(options->cns,options->cni, options->cnd,
							       options->ces,options->cei, options->ced);


  // IPFP used as a refinement method
  IPFPGraphEditDistance<int,int> * algoIPFP = new IPFPGraphEditDistance<int,int>(cf);
  algoIPFP->compute_equiv_mappings(false);
  //algoIPFP->recenterInit();

  // Sinkhorn balanced random init
  if (options->method == string("ipfpe_random_sh")){
    algoIPFP->continuousRandomInit(true);
    options->method = string("ipfpe_multi_random");
  }

  GraphEditDistance<int,int>* ed;
  if(options->method == string("lsape_bunke"))
    ed = new BipartiteGraphEditDistance<int,int>(cf);
  else if( options->method == string("lsape_multi_bunke") )
    ed = new BipartiteGraphEditDistanceMulti<int,int>(cf, options->nep);
  else if( options->method == string("lsape_rw"))
    ed = new RandomWalksGraphEditDistance(cf,options->k);
  else if( options->method == string("lsape_multi_rw") )
    ed = new RandomWalksGraphEditDistanceMulti(cf, options->k, options->nep);
  else if( options->method == string("lsape_multi_greedy") )
    ed = new GreedyGraphEditDistance<int,int>(cf, options->nep);
  //else if( options->method == string("multi_random") )
  //  ed = new RandomInitForIPFP<int,int>(cf, options->nep);

  else if(options->method == string("ipfpe_flat")){
    algoIPFP->continuousFlatInit(true);
    ed = algoIPFP->clone();
  }

  else if(options->method == string("ipfpe_bunke")){
    BipartiteGraphEditDistance<int,int> *ed_init = new BipartiteGraphEditDistance<int,int>(cf);
    ed =new IPFPGraphEditDistance<int,int>(cf,ed_init);
  } else if(options->method == string("ipfpe_multi_bunke")){
    BipartiteGraphEditDistanceMulti<int,int> *ed_init = new BipartiteGraphEditDistanceMulti<int,int>(cf, options->nep);
    ed = new MultistartRefinementGraphEditDistance<int,int>(cf, ed_init, options->nep, algoIPFP);
  // } else if(options->method == string("ipfpe_multi_rw")){
  //   RandomWalksGraphEditDistanceMulti *ed_init = new RandomWalksGraphEditDistanceMulti(cf, options->k, options->nep);
  //   ed = new MultistartRefinementGraphEditDistance<int,int>(cf, ed_init, options->nep, algoIPFP);
  } else if(options->method == string("ipfpe_rw")){
    RandomWalksGraphEditDistance *ed_init = new RandomWalksGraphEditDistance(cf,options->k);
    ed =new IPFPGraphEditDistance<int,int>(cf,ed_init);

  // } else if(options->method == string("ipfpe_multi_random")){
  //   RandomMappingsGED<int,int> *init = new RandomMappingsGED<int,int>();
  //   ed = new MultistartRefinementGraphEditDistance<int,int>(cf, init, options->nep, algoIPFP);

  // } else if(options->method == string("ipfpe_multi_greedy")){
  //   GreedyGraphEditDistance<int,int> *ed_init = new GreedyGraphEditDistance<int,int>(cf, options->nep);
  //   ed = new MultistartRefinementGraphEditDistance<int,int>(cf, ed_init, options->nep, algoIPFP);

  } else if(options->method == string("gnccp")){
    //RandomWalksGraphEditDistance *ed_init = new RandomWalksGraphEditDistance(cf,3 );
    ed = new GNCCPGraphEditDistance<int,int>(cf);//,ed_init);

  } else{
    cerr << "Undefined graph edit distance algorithm "<< endl;
    usage(argv[0]);
    return EXIT_FAILURE;
  }

  ChemicalDataset<double> * dataset = new ChemicalDataset<double>(options->dataset_file.c_str());
  double * distances = computeGraphEditDistance(dataset,ed,options->shuffle, options->nep);

  //Output average distances
  //cout << mean(distances,dataset->size()*dataset->size())<< endl;


  delete algoIPFP;
  delete ed;
  delete dataset;
  delete cf;
  delete options;

  delete [] distances;

  return 0;
}

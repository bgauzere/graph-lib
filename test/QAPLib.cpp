/**
 * @file QAPLib.cpp
 * @author Evariste
 *
 * Programme d'expériences sur QAPLib
 *
 */


// Perform xx time tests to get an average
#define XP_TIME_SAMPLES 1

#include <unistd.h>
#include <iostream>
#include <string>
#include <sstream>

#include <ctime>
#include <chrono>

//#include "xp_output.h"

#include "graph.h"
#include "QAPLibDataset.h"
#include "GraphEditDistance.h"
#include "QAPLibGraph.h"
#include "QAPLibCostFunction.h"
#include "IPFPQAP.h"
#include "MultistartMappingRefinement.h"
#include "RandomMappings.h"

#include "utils.h"

using namespace std;



void usage (char * s)
{
  cerr << "Usage : "<< s << " dataset " << " options"<<endl;
  cerr << "options:" << endl;
  cerr << "\t -m method " << endl;
  cerr << "\t \t Specify the algorithm used to compute edit distance" << endl;
  cerr << "\t -o output_file " << endl;
  cerr << "\t \t Specify a filename for outputing gram matrix" << endl;
  cerr << "\t -c cns,cni,cnd,ces,cei,ced " << endl;
  cerr << "\t \t Specify edit operation costs" << endl;
  cerr << "\t -p n_edit_paths " << endl;
  cerr << "\t \t Specify the number of edit paths to compute GED (lsape_multi)" << endl;
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
  int k = 3;
  int nep = 100; // number of edit paths for allsolution
};

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
    default: /* '?' */
      cerr << "Options parsing failed."  << endl;
      usage(argv[0]);
      exit(EXIT_FAILURE);
    }
  }
  return options;
}




int main (int argc, char** argv)
{

  struct Options * options =   parseOptions(argc,argv);


  options->k = 3;
  QAPLibCost * cf = new QAPLibCost();

  RandomMappings<int,int> *init = new RandomMappings<int,int>();
  IPFPQAP<int,int> * ipfp = new IPFPQAP<int,int>(cf);
  
  MultistartMappingRefinement<int,int> * mipfp = 
        new MultistartMappingRefinement<int,int>(init, options->nep);

  QAPLibDataset * dataset = new QAPLibDataset(options->dataset_file.c_str());
  //double * solutions = computeSolutions<int,int,int>(dataset, ipfp, mipfp, options->shuffle, options->nep);


  if(options->shuffle)
    dataset->shuffleize();

  int N = dataset->size()/2; // The dataset is composed of graph pairs
  double* solutions = new double[N];

  for (int i=0; i<N; i++){
    #ifdef PRINT_TIMES
     clock_t t = clock();
    #endif

    int n = (*dataset)[2*i]->Size();
    int * G1_to_G2 =  new int[n];

    if (options->method == string("ipfp")){
      ipfp->getBetterMapping((*dataset)[2*i], (*dataset)[2*i+1], G1_to_G2, NULL, false);
    }
    else if (options->method == string("mipfp")){
      mipfp->getBestMapping(ipfp, (*dataset)[2*i], (*dataset)[2*i+1], G1_to_G2, NULL);
    }
    
    solutions[i] = ipfp->mappingCost((*dataset)[2*i], (*dataset)[2*i+1], G1_to_G2);


    for (int k=0; k<n; k++){
      cout << G1_to_G2[k]+1 << "  " ;
    } cout << endl;

    delete [] G1_to_G2;

    #ifdef PRINT_TIMES
     t = clock() - t;
     cout << ((float)t) / CLOCKS_PER_SEC << ", " ;
    #endif

    cout << (int)solutions[i];
    cout << endl;

  }
  

  delete ipfp;
  delete mipfp;
  delete init;
  delete dataset;
  delete cf;
  delete options;

  delete [] solutions;

  return 0;
}

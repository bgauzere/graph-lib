/*
 * To compile, run `make` in this folder (/test/).
 */
#include <iostream>
#include <string>
#include <sstream>
#include "../graph_edit_distance.h"

using namespace std;



int main(int argc, char** argv){  

  /* Define a edit distance cost */
  ConstantEditDistanceCost* cf = new ConstantEditDistanceCost(1,3,3,1,3,3);

  /* Load a dataset */
  ChemicalDataset<double>  dataset(argv[1]);

  /* Set some options from the default */
  ged_opts opts = default_refined_opts;
    opts.method = IPFPE_BUNKE;
    opts.ged_output_format = FORMAT_MATRIX;
    opts.dataset_both_dir = true;
    opts.dataset_identity = true;



  /**************************
   * Distances on a dataset *
   **************************/
    
  /* Compute and display distance matrix of the dataset */
  double * matrix = graph_edit_distance<int,int> (dataset, cf, opts);
  for (int i=0; i<dataset.size(); i++){
    for(int j=0; j<dataset.size(); j++){
      cout << matrix[sub2ind(i,j,dataset.size())] << "\t";
    }
    cout << endl;
  }



  /*****************************
   * Mappings between 2 graphs *
   *****************************/
  
  Graph<int,int>& g1 = *(dataset[0]);
  Graph<int,int>& g2 = *(dataset[2]);

  cout << endl;
  cout << " A good mapping from G1 to G2 : " << endl;

  
  /* A best mapping from g1 to g2  */
  ECMapping *mapping = NULL;
  double d = edit_distance_mapping<int,int>( g1, g2, &mapping, NULL, opts);


  /* Display the graph matching
   * Use mapping->f(i) to access assignments of node i of g1 (forward mapping)
   * Use mapping->r(j) to access assignments of node j of graph g2 (reverse mapping)
   * Use mapping->cost() to access the cost of the mapping
   */
  for (uint i=0; i<g1.Size(); i++)
    cout << mapping->f(i) << " ";
  cout << endl;

  for (uint i=0; i<g2.Size(); i++)
    cout << mapping->r(i) << " ";
  cout << endl;

  cout << "Cost : " << d << endl << endl << endl;

  
  
  cout << " 4 mappings from G1 to G2 : " << endl;
  opts.nb_edit_paths = 4;

  
  /* Get a list of refined mappings */
  list<ECMapping*> lmap;
  d = edit_distance_mappings<int,int>( g1, g2, lmap, NULL, opts);

  for (list<ECMapping*>::iterator it = lmap.begin(); it != lmap.end(); it++){
    
    for (uint i=0; i<g1.Size(); i++)
      cout << (*it)->f(i) << " ";
    cout << endl;
    
    for (uint i=0; i<g2.Size(); i++)
      cout <<  (*it)->r(i) << " ";
    cout << endl;
    
    cout << "Cost : " << (*it)->cost()  << endl << endl;
  }


  /* Clean up */
  delete [] matrix;
  delete mapping;
  for (list<ECMapping*>::iterator it = lmap.begin(); it != lmap.end(); it++)
    delete *it;
  
  return EXIT_SUCCESS;
}

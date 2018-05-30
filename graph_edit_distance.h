
#ifndef __GED_FRONTEND_H__
#define __GED_FRONTEND_H__

#include "ECMapping.h"

#include "tinyxml.h"
#include "Dataset.h"
#include "GraphEditDistance.h"
#include "ConstantGraphEditDistance.h"

#include "BipartiteGraphEditDistanceMulti.h"
#include "GreedyGraphEditDistance.h"
#include "RandomWalksGraphEditDistanceMulti.h"
#include "IPFPGraphEditDistance.h"
#include "RandomMappings.h"
#include "MultistartRefinementGraphEditDistance.h"
#include "GNCCPGraphEditDistance.h"

#include "gl_utils.h"


//namespace ged{


  /**
   * @brief Method for graph edit distance approximation.
   */
  enum ged_method {
    LSAPE_BUNKE,    /**< Bipartite based on star assignments cost matrices */
    LSAPE_RW,       /**< Bipartite based on random walks assignments cost matrices */
    LSAPE_GREEDY,   /**< Multi-solution approximating LSAPE_BUNKE */
    IPFPE_BUNKE,    /**< IPFP refining an LSAPE_BUNKE solution */
    IPFPE_RW,       /**< IPFP refining an LSAPE_RW solution */
    IPFPE_GREEDY,   /**< IPFP refining bipartite LSAPE_GREEDY solutions */
    IPFPE_RANDOM,   /**< IPFP with random discrete initializations */
    IPFPE_CONTINUOUS_RANDOM, /**< IPFP with random continuous initializations */
    IPFPE_FLAT,     /**< IPFP with flat continuous initialization */
//  IPFPE_CUSTOM, 
    GNCCP           /**< GCNNP algorithm */
  };

  typedef enum ged_method method_t;


  /**
   * @brief Output format of the graph edit distances computed on a whole dataset
   */
  enum ged_dataset_format {
    FORMAT_DEFAULT,  /**< Default format : matrix if dataset_both_dir option is true */
    FORMAT_VECTOR,   /**< Vectorized verion of the matrix (row first), without uncomputed values */
    FORMAT_MATRIX    /**< Matrix version : `M[sub2ind(i,j,n)]` is the distance from `g_i` to `g_j`*/
  };

  typedef enum ged_dataset_format dataset_format_t;



  /**
   * @brief Options for computing graph edit distance approximations and corresponding mappings
   */
  typedef struct ged_opts{
    method_t method;        /**< Algorithm for approximating the GED */
    int nb_edit_paths;      /**< Number of edit paths for Bipartite multi and IPFP multistart */
    int rw_size;            /**< Size of the random walks (for LSAPE_RW and IPFPE_RW methods) */
    bool dataset_both_dir;  /**< If true, for every pair of graphs \f$\{g_k,g_l\}\f$, compute \f$d(g_k,g_l)\f$ AND \f$d(g_l,g_k)\f$ */
    bool dataset_identity;  /**< If true, for every graph \f$g_i\f$, compute \f$d(g_i,g_i)\f$ */
    bool ipfpe_recenter;    /**< If true, apply a recentering operation in the continuous space to the initialization of IPFP */
    dataset_format_t ged_output_format; /**< Output format of the edit distances of a dataset */
  } opts_t;
  


  /***************************************
   *   Constant default options
   ***************************************/
  
  /**
   * @brief Default options for bipartite assignments
   */
  const opts_t default_bipartite_opts = {
    LSAPE_BUNKE, 32, 3, false, false, false, FORMAT_VECTOR
  };

  /**
   * @brief Default options for refined assignments
   */
  const opts_t default_refined_opts = {
    IPFPE_BUNKE, 16, 3, false, false, false, FORMAT_VECTOR
  };

  /***************************************/



  /**
   * @brief Returns an instance of a graph edit distance method according to the options
   */
  template <class N, class E>
  GraphEditDistance<N,E>*
  ged_from_opts( const opts_t & opts,
		 EditDistanceCost<N,E> * cf )
  {
    IPFPGraphEditDistance<N,E>* algoIPFP = NULL;
    GraphEditDistance<N,E>* ed = NULL;
    BipartiteGraphEditDistanceMulti<N,E>* init_bunke = NULL;
    RandomWalksGraphEditDistanceMulti* init_rw = NULL;
    RandomMappingsGED<N,E>* init_rand = NULL;
    GreedyGraphEditDistance<N,E>* init_greedy = NULL;

    if ( opts.method == IPFPE_BUNKE ||
	 opts.method == IPFPE_RW ||
	 opts.method == IPFPE_GREEDY ||
	 opts.method == IPFPE_FLAT ||
	 opts.method == IPFPE_RANDOM ||
	 opts.method == IPFPE_CONTINUOUS_RANDOM )
      {
	algoIPFP = new IPFPGraphEditDistance<N,E> (cf);
	if (opts.ipfpe_recenter)
	  algoIPFP->recenterInit();
      }


    switch (opts.method) {
    case LSAPE_BUNKE :
      ed = new BipartiteGraphEditDistanceMulti<N,E>(cf, opts.nb_edit_paths);
      break;
    case LSAPE_RW:
      ed = new RandomWalksGraphEditDistanceMulti((ConstantEditDistanceCost*)cf, opts.rw_size, opts.nb_edit_paths);
      break;
    case LSAPE_GREEDY:
      ed = new GreedyGraphEditDistance<N,E>(cf, opts.nb_edit_paths);
      break;

    case IPFPE_FLAT:
      algoIPFP->continuousFlatInit(true);
      ed = algoIPFP->clone();
      break;
      
    case IPFPE_BUNKE:
      init_bunke = new BipartiteGraphEditDistanceMulti<N,E>(cf, opts.nb_edit_paths);
      ed = new MultistartRefinementGraphEditDistance<N,E>(cf, init_bunke, opts.nb_edit_paths, algoIPFP);
      break;
      
    case IPFPE_RW:
      init_rw = new RandomWalksGraphEditDistanceMulti((ConstantEditDistanceCost*)cf,
						      opts.rw_size, opts.nb_edit_paths);
      ed = new MultistartRefinementGraphEditDistance<N,E>(cf, init_rw, opts.nb_edit_paths, algoIPFP);
      break;

    case IPFPE_GREEDY:
      init_greedy = new GreedyGraphEditDistance<N,E>(cf, opts.nb_edit_paths);
      ed = new MultistartRefinementGraphEditDistance<N,E>(cf, init_greedy, opts.nb_edit_paths, algoIPFP);
      break;

    case IPFPE_RANDOM:
      init_rand = new RandomMappingsGED<N,E>();
      ed = new MultistartRefinementGraphEditDistance<N,E>(cf, init_rand, opts.nb_edit_paths, algoIPFP);
      break;

    case IPFPE_CONTINUOUS_RANDOM:
      algoIPFP->continuousRandomInit(true);
      init_rand = new RandomMappingsGED<N,E>();
      ed = new MultistartRefinementGraphEditDistance<N,E>(cf, init_rand, opts.nb_edit_paths, algoIPFP);
      break;

    case GNCCP:
      ed = new GNCCPGraphEditDistance<N,E>(cf);
    }

    return ed;
  }



  /**
   * @brief Compute an approximation of the graph edit distance between `g1` and `g2`
   *
   *  When `cf` is NULL, a constant cost function with the following costs is used :
   * * Node substitution : 1
   * * Node insertion : 3
   * * Node deletion : 3
   * * Edge substitution : 1
   * * Edge insertion : 3
   * * Edge deletion : 3
   *
   * @param  g1 first graph
   * @param  g2 second graph
   * @param  cf cost function. If NULL a constant cost function is used
   * @param  opts options
   */
  template <class N, class E>
  double graph_edit_distance ( Graph<N,E> & g1,
                               Graph<N,E> & g2,
			       EditDistanceCost<N,E> * cf = NULL,
			       const opts_t & opts = default_refined_opts )
  {
    if (!cf){
      cf = new ConstantEditDistanceCost(1,3,3,1,3,3);
    }

    GraphEditDistance<N,E>* ed = ged_from_opts(opts, cf);

    double ged = (*ed)(&g1, &g2);
    delete ed;

    return ged;
  }



  /**
   * @brief Compte an approximation of graph edit distances of graphs in a dataset
   *
   *   The output format corresponds to the one given in the options. When
   *   a matrix is demanded, if `opts.dataset_both_dir` or `opts.dataset_identity` is
   *   false, the corresponding values are uninitialized.
   * 
   *   The output matrix or vector is heap allocated and memory management is left
   *   to the user.
   *
   * @param dataset  A set of graphs
   * @param cf  cost function for the GED. If NULL a constant cost function is used
   * @param opts options
   * @return A matrix or vector of approximations of the GED by the method coresponding to the options
   */
  template <class N, class E, class P>
  double* graph_edit_distance ( Dataset<N,E,P> & dataset,
				EditDistanceCost<N,E> * cf = NULL,
				const opts_t & opts = default_refined_opts )
  {
    if (!cf)
      cf = new ConstantEditDistanceCost(1,3,3,1,3,3);
    GraphEditDistance<N,E>* ed = ged_from_opts(opts, cf);
    
    double* output;
    unsigned int n = dataset.size();
    dataset_format_t format = opts.ged_output_format;

    /* If default output format */
    if (format == FORMAT_DEFAULT && opts.dataset_both_dir) format = FORMAT_MATRIX;
    else if (format == FORMAT_DEFAULT) format = FORMAT_VECTOR;

    /* Allocate output */
    if (format == FORMAT_MATRIX)
      output = new double[n*n];
    else
      if (opts.dataset_identity)
	output = new double[n*(n+1)/2];
      else
	output = new double[n*(n-1)/2];

    /* Compute distances */
    unsigned int j_from_i=1;
    unsigned int index=0;
    if (opts.dataset_both_dir) j_from_i=0;
    for (unsigned int i=0; i<n; i++){
      for (unsigned j=j_from_i*i; j<n; j++){
	if (i==j && !opts.dataset_identity) continue;
	if (format == FORMAT_MATRIX) index = sub2ind(i,j,n);
	else index++;

	output[index] = (*ed)(dataset[i], dataset[j]);
      }
    }

    delete ed;

    return output;
  }


  /** @brief  Compute a graph matching corresponding to the best approximation by a given method
   *
   * @param g1  First graph
   * @param g2  Second graph
   * @param mapping  A pointer to be modified by this function, to store the mapping
   * @param cf  The cost function. If NULL a constant cost function is used
   * @param opts  options
   * @return An approximation of the GED by the method coresponding to the options
   */
  template <class N, class E>
  double edit_distance_mapping ( Graph<N,E> & g1,
				 Graph<N,E> & g2,
				 ECMapping ** mapping,
				 EditDistanceCost<N,E> * cf = NULL,
				 const opts_t & opts = default_refined_opts )
  {
    if (!cf){
      cf = new ConstantEditDistanceCost(1,3,3,1,3,3);
    }

    GraphEditDistance<N,E>* ed = ged_from_opts(opts, cf);

    unsigned int *rho = new unsigned int[g1.Size()];
    unsigned int *varrho = new unsigned int[g2.Size()];

    ed->getOptimalMapping( &g1, &g2, rho, varrho );
    double ged = ed->GedFromMapping(&g1, &g2, rho, g1.Size(), varrho, g2.Size());
    *mapping = new ECMapping(rho, varrho, g1.Size(), g2.Size(), ged);
    
    delete ed;
    return ged;

  }



  template <class N, class E>
  double edit_distance_mappings ( Graph<N,E> & g1,
				  Graph<N,E> & g2,
				  std::list<ECMapping*> & mapping,
				  EditDistanceCost<N,E> * cf = NULL,
				  const opts_t & opts = default_refined_opts )
  {
    /* 1-solution methods : */
    if ( opts.method == GNCCP ){
      ECMapping * map;
      double ged = edit_distance_mapping( g1, g2, &map, cf, opts );
      mapping.push_front( map );
      return ged;
    }

    if (!cf){
      cf = new ConstantEditDistanceCost(1,3,3,1,3,3);
    }

    GraphEditDistance<N,E>* ed = ged_from_opts(opts, cf);
    std::list<unsigned int*> listG1toG2, listG2toG1;

    /* Bipartite-multi or multi-start refinement ? */
    if ( opts.method == LSAPE_BUNKE
	 || opts.method == LSAPE_RW
	 || opts.method == LSAPE_GREEDY ) {
      
      MappingGenerator<int,int>* ed_as_gen = (MappingGenerator<int,int>*) ed;
      std::list<unsigned int*> lsap_maps = ed_as_gen->getMappings(&g1, &g2, opts.nb_edit_paths);
      MultiGed<int,int>::mappings_lsap2lsape(lsap_maps, listG1toG2, listG2toG1, g1.Size(), g2.Size());
    }
    else {
      MultistartRefinementGraphEditDistance<int,int>* ed_as_refinement = dynamic_cast<MultistartRefinementGraphEditDistance<int,int>*>(ed);
      listG1toG2 = ed_as_refinement->getBetterMappings(&g1, &g2);
      listG2toG1 = ed_as_refinement->getReverseMappings();
    }

    /* Store mappings in the output list */
    int minGed = std::numeric_limits<int>::max();
    for ( std::list<unsigned int*>::iterator it = listG1toG2.begin(), jt = listG2toG1.begin();
	  it != listG1toG2.end() && jt != listG2toG1.end();
	  it++, jt++ )
    {
      ECMapping* m = new ECMapping( *it, *jt, g1.Size(), g2.Size(),
				    ed->GedFromMapping(&g1, &g2, *it, g1.Size(), *jt, g2.Size())
				  );
      mapping.push_back(m);
      if (m->cost() < minGed)
	minGed = m->cost();
    }

    delete ed;
    
    return minGed;
  }

//}; // end namespace ged

#endif


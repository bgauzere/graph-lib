# graph-lib

A toolbox for graph kernels and approximate graph edit distance algorithms

## Requirements

#### LSAPE

To use this toolbox, you'll need the [LSAPE library](https://bougleux.users.greyc.fr/lsape/) and set up your 
LSAPE_DIR environment variable. For instance if you use bash, in your .bashrc file :

    export LSAPE_DIR=/path/to/your/local/lsape/include/

or with csh (in your .cshrc):

    setenv LSAPE_DIR /path/to/your/local/lsape/include/

#### TinyXML

This toolbox uses [TinyXML](http://www.grinninglizard.com/tinyxml/) to parse Graphml and GXL graph formats.

#### Eigen

You'll need also the [Eigen library](http://eigen.tuxfamily.org/index.php) in version 3.2 or higher

#### OpenMP

The toolbox uses OpenMP for the multithreaded versions. You can compile the sequential versions without OpenMP.


## How to use the toolbox ?

#### Compile

Requirements : GNU C++ compiler g++ 4.8 or higher.

Open a terminal and move to the proper `graph-lib` directory, then type `make`. This will
create in the directory `bin` the following files :
* A static library `graphlib.a`
* An executable `chemical-edit-distances`

You can use this executable to compute the edit distance of all pair of molecules in a dataset. See the section
Usage below. You can also choose one of the following targets (more options) :

* `with_times` : compute all distances in a dataset and print computation time - **Sequetial version**
* `multithread` : compute all distances in a dataset - **Multithreaded version**
* `multithread_with_times` : compute all distances in a dataset and print computation time - **Multithreaded version**
and type `make <target>` in a terminal, replacing `<target>` by your choice.


#### Usage of the given example

To use the executable `bin/compute-edit-distances` :
     
    ./compute-edit-distances   dataset   [options]  -m  method 

With `dataset` the path to your .ds file listing the graph files of your dataset.
Options can be :
* -s : apply shuffling to the nodes of the graphs
* -p N : number of edit paths set to N (for multiple bipatite and multistart refinement versions)

Methods can be :
(Bipartite)
* **lsape_bunke** - Bipartite based on star assignments cost matrices
* **lsape_multi_bunke** - Multi-solution version of lsape_bunke
* **lsape_rw** - Bipartite based on random walks assignments cost matrices
* **lsape_multi_rw** - Multi-solution version of lsape_rw
* **lsape_multi_greedy** - Multi-solution approximating lsape_bunke

(Refinements)
* **ipfpe_flat** - IPFP with flat continuous initialization
* **ipfpe_bunke** - IPFP refining an lsape_bunke solution
* **ipfpe_multi_bunke** - Multistart IPFP refining bipartite lsape_multi_bunke solutions
* **ipfpe_rw** - IPFP refining an lsape_rw solutions
* **ipfpe_multi_rw** - Multistart IPFP refining bipartite lsape_multi_rw solutions
* **ipfpe_multi_greedy** - Multistart IPFP refining bipartite lsape_multi_greedy solutions
* **ipfpe_multi_random** - Multistart IPFP with random discrete initializations
* **gnccp** - GNCCP algorithm


#### Use the toolbox in your code

If you want to use this toolbox in your code, you can include the file `graph_edit_distance.h` in the root
directory. This header defines 4 helper functions to compute approximate graph edit distances and graph matchings :

* `double  graph_edit_distance ( Graph&,  Graph&,  EditDistanceCost*,  options )`
* `double* graph_edit_distance ( Dataset&  EditDistanceCost*,  options )`
* `double  edit_distance_mapping ( Graph&,  Graph&,  ECMapping**,  EditDistanceCost*,  options )`
* `double  edit_distance_mappings ( Graph&,  Graph&,  std::list<ECMapping*>,  EditDistanceCost*,  options )`

The first one and the second one allow to compute approximate graph edit distance between resp. two graphs and
all pairs of a dataset. The last ones compute graph mappings. Check out the documentation in this file, as well as 
the example file : `test/frontend_example.cpp`.

Basically you will need to create at least two graphs with the `Graph<NodeAttr, EdgeAttr>` class, or load a graph
dataset with the `Dataset<NodeAttr, EdgeAttr, PropertyType>` class, a cost function adapted to your graphs (take a
look at `include/GraphEditDistance.h` for more details),  and set up some options.
Then you can call one of the later functions :

```
  /* Define a edit distance cost */
  ConstantEditDistanceCost* cf = new ConstantEditDistanceCost(1,3,3,1,3,3);

  /* Load a dataset */
  ChemicalDataset<double>  dataset(argv[1]);

  /* Set some options from the default */
  ged_opts opts = default_refined_opts;
  opts.method = IPFPE_RW;
  opts.ged_output_format = FORMAT_MATRIX;
  opts.dataset_both_dir = true;
  opts.dataset_identity = true;
  
  double * matrix = graph_edit_distance<int,int> (dataset, cf, opts);
```

The mappings are encoded as two arrays encapsulated in the class `ECMapping`. The first array corresponds to the 
mapping of the nodes of g1 and is accessible via `ECMapping::f(unsigned int)`. The second one corresponds to the mapping of 
the nodes of g2 and is accessible via `ECMapping::r(unsigned int)`. The cost of the mapping is accessible via 
`ECMapping::cost()`

```
  /* A best mapping from g1 to g2  */
  ECMapping *mapping = NULL;
  double d = edit_distance_mapping<int,int>( g1, g2, &mapping, NULL, opts);
  
  for (uint i=0; i<g1.Size(); i++) /* Show the mapping of g1's nodes */
    cout << mapping->f(i) << " ";
  cout << endl;
```



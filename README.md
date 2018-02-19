# graph-lib

A toolbox for graph kernels and approximate graph edit distance algorithms

## Requirements

#### LSAPE

To use this toolbox, you'll need the [LSAPE library](https://bougleux.users.greyc.fr/lsape/) and set up your 
LSAPE_DIR environment variable. For instance if you use bash, in your .bashrc file :

    export LSAPE_DIR /path/to/your/local/lsape/include/

or with csh (in your .cshrc):

    setenv LSAPE_DIR /path/to/your/local/lsape/include/

#### TinyXML

This toolbox uses [TinyXML](http://www.grinninglizard.com/tinyxml/) to parse Graphml and GXL graph formats.

#### Eigen

You'll need also the [Eigen library](http://eigen.tuxfamily.org/index.php) in version 3.2 or higher

#### OpenMP

The toolbox uses OpenMP for the multithreaded versions. You can compile the sequential versions without OpenMP.

## Compile

Requirements : GNU C++ compiler g++ 4.8 or higher.

Choose one of the following rules given what you want to compile :
* with_times : compute all distances in a dataset and print computation time - **Sequetial version**
* multithread : compute all distances in a dataset - **Multithreaded version**
* multithread_with_times : compute all distances in a dataset and print computation time - **Multithreaded version**
and type `make <rule>` in a terminal, replacing <rule> by your choice.


## Usage

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

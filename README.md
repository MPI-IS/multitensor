MultiTensor
===========

The **MutliTensor** is a library for multilayer network tensor factorization that can be used
for community detection, link prediction and measure layer interdependence.


Installation
============

The project uses [CMake](https://cmake.org/) for building:

- the main library
- the [Python](https://www.python.org/) bindings
- the tests
- the documentation


Requirements
------------

- [Boost](http://www.boost.org) for the unit tests
- [CMake](https://cmake.org/) for building the project


Configuring and building the project with CMake
-----------------------------------------------

To configure and build the project, simply type in:

```
cd $MULTI_TENSOR_DIR
mkdir build
cd build
cmake -DBOOST_ROOT=$BOOST_PATH ..
make
```
where `Boost` is installed under `$BOOST_PATH` and `$MULTI_TENSOR_DIR` is the entry point of the library repository.


Python bindings (optional)
--------------------------
Building the Python bindings requires `numpy`. Although it is not mandatory, we suggest to also
use `cython` to regenerate the library source file.
These packages can be installed using [pip](https://pypi.python.org/pypi/pip):

```
pip install numpy cython
```

Then, run your `CMake` with the option `ENABLE_PYTHON_WRAPPER` set to **ON**:

```
cmake -DBOOST_ROOT=$BOOST_PATH -ENABLE_PYTHON_WRAPPER=ON ..
```

**Note**: you should use a Python virtual environment, and `CMake` should be run in this virtual environment.


Running the tests
-----------------
Once the project is built, just type

```
make test
```

This will run the unit tests, also for Python if the extension is enabled.


Usage
=====

## What's included:
- `cpp` : c++ version of the algorithm. Faster than the Python one.
- `python` : Python version. Slower than c++.
- `data` : Contains sample adjacency files to test the code.

## Requirements:
Need to make a directory called `data` outside the `cpp` and `python` folders. 
To make one, just type from the command line, inside that folder: 
* `mkdir data`

## Input format.
The multilayer adjacency matrix should be formatted as an edge list with L+3 columns:

`E node1 node2 3 0 0 1`

The first columns tells the algorithm that the row denotes an edge; the second and third are the source and target nodes of that edge, respectively; l+3 column tells if there is that edge in the l-th layer and the weigth (must be integer). In this example the edge node1 --> node2 exists in layer 1 with weight 3 and in layer 4 with weight 1, but not in layer 2 and 3.

Note: if the network is undirected, you only need to input each edge once. You then need to specificy to the algotihm that you are considering the undirected case: for the `cpp` version this is done by running `./MultiTensor_undirected` (first you need to compile it by changing the Makefile accordingly); for the `python` version this is done by giving as a command line input parameter `-u=1`. 

## Output.
Three files will be generated inside the `data` folder: the two NxK membership matrices `U` and `V`, and the KxK layer affinity matrix `W`. Supposing that K=4 and `E=".dat"` the output files will be inside `data` folder with names:
- `u_K4.dat`
- `v_K4.dat`
- `w_K4.dat`

The first line outputs the Max Likelihood among the realizations.
For the membership files, the subsequent lines contain L+1 columns: the first one is the node label, the follwing ones are the (not normalized) membership vectors' entries.
For the affinity matrix file, the subsequent lines start with the number of the layer and then the matrix for that layer.
For the restricted assortative version only the diagonal entries of the affinity matrix are printed. The first entry of each row is the layer index.


References
==========

The **MutliTensor** implements the algorithm described in:

De Bacco, C., Power, E. A., Larremore, D. B., & Moore, C. (2017). *Community detection, link prediction, and layer interdependence in multilayer networks.* Physical Review E, 95(4), 042317.

If you use this code please cite this [article](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.95.042317).  
A _preprint_ version can be found [here](http://cdebacco.com/files/multitensor.pdf) or [here](https://arxiv.org/abs/1701.01369).


Authors
=======

[Caterina De Bacco](caterina.debacco@tuebingen.mpg.de)


[Jean-Claude Passy](jean-claude.passy@tuebignen.mpg.de)


License
=======

GNU GPL version 3 (see LICENSE.md)


Copyright
=========

(c) 2019, Max Planck Society / Software Workshop - Max Planck Institute for Intelligent Systems



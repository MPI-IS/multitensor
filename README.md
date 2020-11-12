MultiTensor
===========

The **MutliTensor** is a library for multilayer network tensor factorization that can be used
for community detection, link prediction, and to measure layer interdependence.


Installation
============

The project is written in `C++` and includes

- the main library
- the Python extension
- the tests
- the benchmark
- the documentations


Requirements
------------

- [Boost libraries](http://www.boost.org)
- [CMake](https://cmake.org/) for building the project
- [Python 3.6+](https://www.python.org/) and [Numpy](https://numpy.org/)
for building the Python bindings and running the benchmark
- [Doxygen](https://www.doxygen.nl/index.html) for generating the main documentation


Configuring and building the project with CMake
-----------------------------------------------

To configure and build the project, simply type in:
```
$ cd $MULTI_TENSOR_DIR
$ mkdir build
$ cd build
$ cmake -DBOOST_ROOT=$BOOST_PATH -DCMAKE_BUILD_TYPE=$BUILD_TYPE ..
$ make
```
where
* `$MULTI_TENSOR_DIR` is the root of this repository
* `$BOOST_PATH` is the location of the `Boost` libraries (e.g. `/usr/local/Cellar/boost/1.72.0_3/`)
* `$BUILD_TYPE` is the configuration you wish to build (usually `Release` or `Debug`)

Documentation
-------------

The main documentation is written with Doxygen. It can be built the following way:
```
$ make doxygen
```

Running the tests
-----------------
Once the project is built, just type
```
$ make test
```

This will run the unit tests and the functional tests,
also for Python if the extensions are enabled.

Running the benchmark
---------------------

To run the benchmark (Table III) from the [De Bacco et al. 2017](#References),
run from the `build` directory:
```
$ cd benchmark
$ python benchmark_runner.py --help
$ python benchmark_runner.py --test $TEST_NAME
```

Usage
=====

C++ API
-------

A binary `Multitensor` is provided to use the C++ implementation.
Use the `help` option to obtain more details about the binary:
```
$ ./Multitensor --help
```

You can run it against the examples provided in the `data` folder
in the source repository, for instance:
```
$ ./Multitensor --a $MULTITENSOR_SRC_DIR/data/main/adjacency.dat --k 2
```
By default, the results of the algorithm are written in the `results` folder.

Input format
------------

The multilayer adjacency matrix should be formatted as an edge list with L+3 columns:

```
E node1 node2 3 0 0 1
```

* The first column (`E`) tells the algorithm that the row denotes an edge
* The second and third are the labels source and target nodes of that edge, respectively
* The remaining columns represent the edge weight in each layer

In the example above the edge *node1* --> *node2* exists in layer 1 with weight 3
and in layer 4 with weight 1, but not in layer 2 and 3.

Output format
-------------

Four files will be generated inside an output folder (by default named `results`):`
* some information about the run of the algorithm
* the NxK membership matrix `U`
* the NxK membership matrix `V`
* the KxK layer affinity matrix `W`

For the membership matrices, the first column is the node label and the following ones the membership vector entries.


As for the affinity matrix, it is organized in blocks separated by an empty line.
Each block starts with the layer number followed by the matrix for that layer. For the assortative version only the diagonal entries of the affinity matrix are printed. The first entry of each row is the layer index.

Python API
----------

A Python extension build with [Cython](https://cython.org/) is also provided.
To build it:
```
$ pip install -U numpy cython # cython is optional
$ make multitensor_py
```

The API is described in the Sphinx documentation
which can be built the following way:
```
$ pip install -U sphinx sphinx-bootstrap-theme
$ make sphinx
```


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

[Ivan Oreshnikov](ivan.oreshnikov@tuebignen.mpg.de)

License
=======

GNU GPL version 3 (see LICENSE.md)


Copyright
=========

(c) 2019, Max Planck Society / Software Workshop - Max Planck Institute for Intelligent Systems

# This is an auxiliary module providing the reader functions. Since
# the input format might not be the standard in the research
# community, there is a chance that this module will not survive until
# the public release. Meanwhile we use it in to load the test data.


from libcpp.vector cimport vector
from libcpp.string cimport string

import numpy
cimport numpy


cdef extern from "boost/filesystem.hpp" namespace "boost::filesystem":
    cppclass path:
        pass


cdef extern from "app_utils.hpp":
    cdef void read_adjacency_data[weight_t](
        const path &filename,
        vector[size_t] &edges_start,
        vector[size_t] &edges_end,
        vector[weight_t] &edges_weights
    )

    cdef void read_affinity_data(
        const path &filename,
        const bint &assortative,
        vector[numpy.float_t] &w
    )


def py_read_adjacency_data(string filename):
    # We define empty vectors that will be populated by the vertices
    # correponding to the beginning and end the end of the network
    # edges, as well as a vector for the unraveled edge weight matrix.
    cdef vector[size_t] edges_start
    cdef vector[size_t] edges_end
    cdef vector[numpy.float_t] edges_weights

    # We read the data using the C function.
    read_adjacency_data[numpy.float_t](
        <const path &> filename,
        <vector[size_t] &> edges_start,
        <vector[size_t] &> edges_end,
        <vector[numpy.float_t] &> edges_weights,
    )

    # And then we copy it into the numpy arrays.
    py_edges_start = numpy.zeros([edges_start.size()], dtype=numpy.int)
    py_edges_end = numpy.zeros([edges_start.size()], dtype=numpy.int)
    py_edges_weights = numpy.zeros([edges_weights.size()], dtype=numpy.float)

    for i in range(edges_start.size()):
        py_edges_start[i] = edges_start[i]

    for i in range(edges_end.size()):
        py_edges_end[i] = edges_end[i]

    for i in range(edges_weights.size()):
        py_edges_weights[i] = edges_weights[i]

    return py_edges_start, py_edges_end, py_edges_weights


def py_read_affinity_data(
    string filename,
    size_t nof_groups,
    size_t nof_layers,
    bint assortative
):
    # We define the vector for the affinity matrix.
    cdef vector[numpy.float_t] w

    # And resize it to the appropriate size depending on whether we're
    # working with the assortative case or not.
    if assortative:
        w.resize(nof_groups * nof_layers)
    else:
        w.resize(nof_groups * nof_groups * nof_layers)

    # Read the data using the C fucntion.
    read_affinity_data(
        <const path &> filename,
        <const bint &> assortative,
        <vector[numpy.float_t] &> w
    )

    # Copy the result into the numpy array.
    py_w = numpy.zeros([w.size()], dtype=numpy.float)
    for i in range(w.size()):
        py_w[i] = w[i]

    return py_w

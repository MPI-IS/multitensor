# Copyright (c) 2019, Max Planck Society / Software Workshop - Max Planck Institute for Intelligent Systems
# Distributed under the GNU GPL license version 3
# See file LICENSE.md or at https://github.com/MPI-IS/multitensor/LICENSE.md

'''
Python wrapper for the multitensor library

:author: Ivan Oreshnikov <ivan.oreshnikov@tuebingen.mog.de>
:author: Jean-Claude Passy <jean-claude.passy@tuebingen.mpg.de>
'''

from libcpp.string cimport string
from libc.time cimport time, time_t
from libcpp.vector cimport vector

# To dereference pointers
# http://cython.readthedocs.io/en/latest/src/userguide/wrapping_CPlusPlus.html#c-operators-not-compatible-with-python-syntax
from cython.operator cimport dereference as deref

import logging
import numpy
cimport numpy

# 0 - Imports
cdef extern from "multitensor/utils.hpp" namespace "multitensor::utils":
    cdef cppclass Report:
        size_t nof_realizations
        vector[double] vec_L2
        vector[size_t] vec_iter
        vector[const char *] vec_term_reason
        double duration

        Report() except +
        double max_L2()

cdef extern from "multitensor/tensor.hpp" namespace "multitensor::tensor":
    cdef cppclass Tensor[scalar_t]:
        size_t nrows
        size_t ncols
        size_t ntubes
        vector[scalar_t] data

        Tensor() except +
        Tensor(size_t nrows, size_t ncols, size_t ntubes=1) except +

        void resize(size_t nrows_, size_t ncols_, size_t ntubes)
        size_t get_nrows()
        size_t get_ncols()
        size_t get_ntubes()

        scalar_t & operator()(const size_t i, const size_t j, const size_t k)

    cdef cppclass Matrix[scalar_t](Tensor[scalar_t]):
        Matrix() except +
        Matrix(size_t nrows, size_t ncols) except +

        void resize(size_t nrows, size_t ncols)
        scalar_t & operator()(const size_t i, const size_t j)

    cdef cppclass SymmetricTensor[scalar_t](Tensor[scalar_t]):
        SymmetricTensor() except +
        SymmetricTensor(size_t nrows, size_t ntubes) except +

        void resize(size_t nrows, size_t ntubes)

    cdef cppclass DiagonalTensor[scalar_t](Tensor[scalar_t]):
        DiagonalTensor() except +
        DiagonalTensor(size_t nrows, size_t ntubes)


# Algorithm
# Here we define all the necessary types.
# We start with exporting direction selectors from boost graph library.
cdef extern from "boost/graph/graph_selectors.hpp" namespace "boost":
    struct directedS:
        pass
    struct undirectedS:
        pass
    struct bidirectionalS:
        pass

# We are interested in only two cases -- either undirected or
# bidirectional graphs. The concrete implementation of the algorithm
# is chosen by a template parameter, hence the fused type.
ctypedef fused direction_t:
    undirectedS
    bidirectionalS

# We might end up using either diagonal or symmetric tensors for the
# affinity data. This is a template parameter of the algorithm.
ctypedef fused affinity_t:
    DiagonalTensor[numpy.float_t]
    SymmetricTensor[numpy.float_t]

# We have two separate modes of affinity data initialization -- either
# as a random tensor or from the inital data loaded from disk. This is
# a template parameter of the algorithm.
cdef extern from "multitensor/initialization.hpp" namespace "multitensor::initialization":
    cppclass init_symmetric_tensor_random:
        pass
    cppclass init_symmetric_tensor_from_initial[affinity_t]:
        pass

# At the moment we support only the integer numbers as the graph
# vertices. We might want to extend this later.
ctypedef numpy.int_t vertex_t

# Edge weights can be defined either by an integer number or by a
# floating point number. This is a template parameter of the
# algorithm.
ctypedef fused weight_t:
    numpy.int_t
    numpy.float_t


# This is a python wrapper for the report returned by the solver.
cdef class ReportWrapper:
    """Wrapper for the report returned by the solver."""

    cdef Report c_obj

    @property
    def nof_realizations(self):
        """Number of realizations."""
        return self.c_obj.nof_realizations

    @nof_realizations.setter
    def nof_realizations(self, nof_realizations):
        self.c_obj.nof_realizations = nof_realizations

    @property
    def vec_L2(self):
        """Likelihood for each realization."""
        return self.c_obj.vec_L2

    @vec_L2.setter
    def vec_L2(self, vec_L2):
        self.c_obj.vec_L2 = vec_L2

    @property
    def vec_iter(self):
        """Number of iterations for each realization."""
        return self.c_obj.vec_iter

    @vec_iter.setter
    def vec_iter(self, vec_iter):
        self.c_obj.vec_iter = vec_iter

    @property
    def vec_term_reason(self):
        """Reason for terminating the solver for each realization."""
        return self.c_obj.vec_term_reason

    @vec_term_reason.setter
    def vec_term_reason(self, vec_term_reason):
        self.c_obj.vec_term_reason = vec_term_reason

    @property
    def duration(self):
        """Duration (in seconds) of the full run."""
        return self.c_obj.duration

    @duration.setter
    def duration(self, duration):
        self.c_obj.duration = duration

    @property
    def max_L2(self):
        """Maximum likelihood."""
        return self.c_obj.max_L2()


# A small utility function for calculating the number of network
# vertices from the edge endpoints.
cdef extern from "multitensor/utils.hpp" namespace "multitensor::utils":
    cdef size_t get_num_vertices[vertex_t](
        const vector[vertex_t] & edges_start,
        const vector[vertex_t] & edges_end
    )

    cdef cppclass RandomGenerator[R, D]:
        RandomGenerator(time_t seed) except +


cdef extern from "<random>" namespace "std":
    cdef cppclass mt19937_64 "std::mt19937_64":
        pass
    cdef cppclass uniform_real_distribution "std::uniform_real_distribution<double>":
        pass


# Main entry point of the algorithm.
cdef extern from "multitensor/main.hpp" namespace "multitensor":
    cdef Report c_multitensor_factorization "multitensor::multitensor_factorization"[
        direction_t, affinity_t, affinity_init_t, vertex_t, weight_t](
        const vector[vertex_t] & edges_start,
        const vector[vertex_t] & edges_end,
        const vector[weight_t] & edges_weight,
        const size_t & nof_realizations,
        const size_t & max_nof_iterations,
        const size_t & nof_convergences,
        vector[vertex_t] & labels,
        Matrix[numpy.float_t] & u,
        Matrix[numpy.float_t] & v,
        vector[numpy.float_t] & affinity
    ) except +RuntimeError


def run(
    adjacency_filename,
    nof_groups,
    directed=True,
    assortative=False,
    nof_realizations=1,
    max_nof_iterations=500,
    nof_convergences=10,
    init_affinity_filename=None,
    weigths_dtype=float,
    seed=None
):
    """
    Runs the multitensor factorization algorithm.

    :param str adjacency_filename: Name of the file containing the adjacency data
    :param int nof_groups: Number of groups
    :param bool directed: Whether the network is directed (True) or undirected (False)
    :param bool assortative: If True, assumes an assortative model
    :param int nof_realizations: Number of realizations
    :param int max_nof_iterations: Maximum number of iterations in each realization
    :param int nof_convergences: Number of successive passed convergence criteria
        for declaring the results converged
    :param str init_affinity_filename: Name of the file containing the initial affinity data
    :param dtype weigths_dtype: Type used for the edge weights
    :param int seed: Seed for the random generator (mt19937 with uniform distribution)
    :returns: 4-tuple containing:

        * the numpy array linking outgoing vertices
        * the numpy array linking incoming vertices
        * the numpy array containing the affinity values
        * the detailed report
    :rtype: tuple(numpy.array, numpy.array, numpy.array, ReportWrapper)
    """
    # Load adjacency file
    adj_data = numpy.loadtxt(adjacency_filename)

    edges_start = adj_data[:, 0].astype(int)
    edges_end = adj_data[:, 2].astype(int)
    edges_weights = adj_data[:, 2:].astype(weigths_dtype).ravel()

    # The underlying C function deduces the number of groups from the
    # shane of the affinity matrix passed as an input. For a python
    # version we want to keep a more python approach and pass an
    # explicit parameter. This implies that we have to repeat the math
    # done in C function in reverse.
    cdef size_t nof_edges = edges_start.size
    cdef size_t nof_vertices = get_num_vertices[vertex_t](
        < const vector[vertex_t] & > edges_start,
        < const vector[vertex_t] & > edges_end)
    cdef size_t nof_layers = edges_weights.size / nof_edges

    # Preallocate the output matrices.
    cdef Matrix[numpy.float_t] c_u = Matrix[numpy.float_t](nof_vertices, nof_groups)
    cdef Matrix[numpy.float_t] c_v = Matrix[numpy.float_t](0, 0)

    # Preallocate the affinity vector depending on the
    # assortative/nonassortative case.
    cdef vector[numpy.float_t] c_affinity

    # Initialize the affinity vector
    if init_affinity_filename:
        # read the affinity file
        w_data = numpy.loadtxt(init_affinity_filename)

        if assortative:
            init_affinity = w_data[:, 1:].ravel()
        else:
            init_affinity = (numpy.diag(l) for l in w_data[:, 1:])
            init_affinity = numpy.concatenate([l.ravel() for l in init_affinity])
        c_affinity = < vector[numpy.float_t] > init_affinity
    else:
        if assortative:
            affinity_size = nof_groups * nof_layers
        else:
            affinity_size = nof_groups * nof_groups * nof_layers
        c_affinity = vector[numpy.float_t](< size_t > affinity_size)

    # Create an empty label vector. Not sure why we need this :)
    cdef vector[vertex_t] labels = vector[vertex_t](nof_vertices)

    # Detailed report
    report = ReportWrapper()

    # Random generator
    seed = seed if seed is not None else time(NULL)

    # Cannot stack-allocate C++ objects with constructor arguments in cython
    # See https://stackoverflow.com/questions/45991342/directly-call-c-struct-constructor-from-cython
    # and http://cython.readthedocs.io/en/latest/src/userguide/wrapping_CPlusPlus.html
    # So here we use a pointer
    cdef RandomGenerator[mt19937_64, uniform_real_distribution] * rng = \
        new RandomGenerator[mt19937_64, uniform_real_distribution](seed)

    # The code below reproduces the switch case from MultiTensor.hpp
    # (starting on line 175 as of moment of writing). This is done for
    # the case of weight_t being an int.
    try:
        if weigths_dtype is int and not directed and not assortative and not init_affinity_filename:
            # case 0: undirected + non-assortative + w random
            report.c_obj = c_multitensor_factorization[
                undirectedS,
                SymmetricTensor[numpy.float_t],
                init_symmetric_tensor_random,
                vertex_t,
                numpy.int_t
            ](
                < const vector[vertex_t] & > edges_start,
                < const vector[vertex_t] & > edges_end,
                < const vector[numpy.int_t] & > edges_weights,
                nof_realizations, max_nof_iterations, nof_convergences,
                labels, c_u, c_v, c_affinity, deref(rng)
            )

        if weigths_dtype is int and directed and not assortative and not init_affinity_filename:
            # case 1: directed + non-assortative + w random
            c_v.resize(nof_vertices, nof_groups)
            report.c_obj = c_multitensor_factorization[
                bidirectionalS,
                SymmetricTensor[numpy.float_t],
                init_symmetric_tensor_random,
                vertex_t,
                numpy.int_t
            ](
                < const vector[vertex_t] & > edges_start,
                < const vector[vertex_t] & > edges_end,
                < const vector[numpy.int_t] & > edges_weights,
                nof_realizations, max_nof_iterations, nof_convergences,
                labels, c_u, c_v, c_affinity, deref(rng)
            )

        if weigths_dtype is int and not directed and assortative and not init_affinity_filename:
            # case 2: undirected + assortative + w random
            report.c_obj = c_multitensor_factorization[
                undirectedS,
                DiagonalTensor[numpy.float_t],
                init_symmetric_tensor_random,
                vertex_t,
                numpy.int_t
            ](
                < const vector[vertex_t] & > edges_start,
                < const vector[vertex_t] & > edges_end,
                < const vector[numpy.int_t] & > edges_weights,
                nof_realizations, max_nof_iterations, nof_convergences,
                labels, c_u, c_v, c_affinity, deref(rng)
            )

        if weigths_dtype is int and directed and assortative and not init_affinity_filename:
            # case 3: directed + assortative + w random
            c_v.resize(nof_vertices, nof_groups)
            report.c_obj = c_multitensor_factorization[
                bidirectionalS,
                DiagonalTensor[numpy.float_t],
                init_symmetric_tensor_random,
                vertex_t,
                numpy.int_t
            ](
                < const vector[vertex_t] & > edges_start,
                < const vector[vertex_t] & > edges_end,
                < const vector[numpy.int_t] & > edges_weights,
                nof_realizations, max_nof_iterations, nof_convergences,
                labels, c_u, c_v, c_affinity, deref(rng)
            )

        if weigths_dtype is int and not directed and not assortative and init_affinity_filename:
            # case 4: undirected + non-assortative + w from file
            report.c_obj = c_multitensor_factorization[
                undirectedS,
                SymmetricTensor[numpy.float_t],
                init_symmetric_tensor_from_initial[SymmetricTensor[numpy.float_t]],
                vertex_t,
                numpy.int_t
            ](
                < const vector[vertex_t] & > edges_start,
                < const vector[vertex_t] & > edges_end,
                < const vector[numpy.int_t] & > edges_weights,
                nof_realizations, max_nof_iterations, nof_convergences,
                labels, c_u, c_v, c_affinity, deref(rng)
            )

        if weigths_dtype is int and directed and not assortative and init_affinity_filename:
            # case 5: directed + non-assortative + w from file
            c_v.resize(nof_vertices, nof_groups)
            report.c_obj = c_multitensor_factorization[
                bidirectionalS,
                SymmetricTensor[numpy.float_t],
                init_symmetric_tensor_from_initial[SymmetricTensor[numpy.float_t]],
                vertex_t,
                numpy.int_t
            ](
                < const vector[vertex_t] & > edges_start,
                < const vector[vertex_t] & > edges_end,
                < const vector[numpy.int_t] & > edges_weights,
                nof_realizations, max_nof_iterations, nof_convergences,
                labels, c_u, c_v, c_affinity, deref(rng)
            )

        if weigths_dtype is int and not directed and assortative and init_affinity_filename:
            # case 6: undirected + assortative + w from file
            report.c_obj = c_multitensor_factorization[
                undirectedS,
                DiagonalTensor[numpy.float_t],
                init_symmetric_tensor_from_initial[DiagonalTensor[numpy.float_t]],
                vertex_t,
                numpy.int_t
            ](
                < const vector[vertex_t] & > edges_start,
                < const vector[vertex_t] & > edges_end,
                < const vector[numpy.int_t] & > edges_weights,
                nof_realizations, max_nof_iterations, nof_convergences,
                labels, c_u, c_v, c_affinity, deref(rng)
            )

        if weigths_dtype is int and directed and assortative and init_affinity_filename:
            # case 7: directed + assortative + w from file
            c_v.resize(nof_vertices, nof_groups)
            report.c_obj = c_multitensor_factorization[
                bidirectionalS,
                DiagonalTensor[numpy.float_t],
                init_symmetric_tensor_from_initial[DiagonalTensor[numpy.float_t]],
                vertex_t,
                numpy.int_t
            ](
                < const vector[vertex_t] & > edges_start,
                < const vector[vertex_t] & > edges_end,
                < const vector[numpy.int_t] & > edges_weights,
                nof_realizations, max_nof_iterations, nof_convergences,
                labels, c_u, c_v, c_affinity, deref(rng)
            )

        # At the moment Cython does not support using fused types to
        # instantiate a template argument, and we have to explicitly pick
        # the implementation for every type variant. This is the same
        # switch case repeated for floating-piont weight_t.
        if weigths_dtype is float and not directed and not assortative and not init_affinity_filename:
            # case 0: undirected + non-assortative + w random
            report.c_obj = c_multitensor_factorization[
                undirectedS,
                SymmetricTensor[numpy.float_t],
                init_symmetric_tensor_random,
                vertex_t,
                numpy.float_t
            ](
                < const vector[vertex_t] & > edges_start,
                < const vector[vertex_t] & > edges_end,
                < const vector[numpy.float_t] & > edges_weights,
                nof_realizations, max_nof_iterations, nof_convergences,
                labels, c_u, c_v, c_affinity, deref(rng)
            )

        if weigths_dtype is float and directed and not assortative and not init_affinity_filename:
            # case 1: directed + non-assortative + w random
            c_v.resize(nof_vertices, nof_groups)
            report.c_obj = c_multitensor_factorization[
                bidirectionalS,
                SymmetricTensor[numpy.float_t],
                init_symmetric_tensor_random,
                vertex_t,
                numpy.float_t
            ](
                < const vector[vertex_t] & > edges_start,
                < const vector[vertex_t] & > edges_end,
                < const vector[numpy.float_t] & > edges_weights,
                nof_realizations, max_nof_iterations, nof_convergences,
                labels, c_u, c_v, c_affinity, deref(rng)
            )

        if weigths_dtype is float and not directed and assortative and not init_affinity_filename:
            # case 2: undirected + assortative + w random
            report.c_obj = c_multitensor_factorization[
                undirectedS,
                DiagonalTensor[numpy.float_t],
                init_symmetric_tensor_random,
                vertex_t,
                numpy.float_t
            ](
                < const vector[vertex_t] & > edges_start,
                < const vector[vertex_t] & > edges_end,
                < const vector[numpy.float_t] & > edges_weights,
                nof_realizations, max_nof_iterations, nof_convergences,
                labels, c_u, c_v, c_affinity, deref(rng)
            )

        if weigths_dtype is float and directed and assortative and not init_affinity_filename:
            # case 3: directed + assortative + w random
            c_v.resize(nof_vertices, nof_groups)
            report.c_obj = c_multitensor_factorization[
                bidirectionalS,
                DiagonalTensor[numpy.float_t],
                init_symmetric_tensor_random,
                vertex_t,
                numpy.float_t
            ](
                < const vector[vertex_t] & > edges_start,
                < const vector[vertex_t] & > edges_end,
                < const vector[numpy.float_t] & > edges_weights,
                nof_realizations, max_nof_iterations, nof_convergences,
                labels, c_u, c_v, c_affinity, deref(rng)
            )

        if weigths_dtype is float and not directed and not assortative and init_affinity_filename:
            # case 4: undirected + non-assortative + w from file
            report.c_obj = c_multitensor_factorization[
                undirectedS,
                SymmetricTensor[numpy.float_t],
                init_symmetric_tensor_from_initial[SymmetricTensor[numpy.float_t]],
                vertex_t,
                numpy.float_t
            ](
                < const vector[vertex_t] & > edges_start,
                < const vector[vertex_t] & > edges_end,
                < const vector[numpy.float_t] & > edges_weights,
                nof_realizations, max_nof_iterations, nof_convergences,
                labels, c_u, c_v, c_affinity, deref(rng)
            )

        if weigths_dtype is float and directed and not assortative and init_affinity_filename:
            # case 5: directed + non-assortative + w from file
            c_v.resize(nof_vertices, nof_groups)
            report.c_obj = c_multitensor_factorization[
                bidirectionalS,
                SymmetricTensor[numpy.float_t],
                init_symmetric_tensor_from_initial[SymmetricTensor[numpy.float_t]],
                vertex_t,
                numpy.float_t
            ](
                < const vector[vertex_t] & > edges_start,
                < const vector[vertex_t] & > edges_end,
                < const vector[numpy.float_t] & > edges_weights,
                nof_realizations, max_nof_iterations, nof_convergences,
                labels, c_u, c_v, c_affinity, deref(rng)
            )

        if weigths_dtype is float and not directed and assortative and init_affinity_filename:
            # case 6: undirected + assortative + w from file
            report.c_obj = c_multitensor_factorization[
                undirectedS,
                DiagonalTensor[numpy.float_t],
                init_symmetric_tensor_from_initial[DiagonalTensor[numpy.float_t]],
                vertex_t,
                numpy.float_t
            ](
                < const vector[vertex_t] & > edges_start,
                < const vector[vertex_t] & > edges_end,
                < const vector[numpy.float_t] & > edges_weights,
                nof_realizations, max_nof_iterations, nof_convergences,
                labels, c_u, c_v, c_affinity, deref(rng)
            )

        if weigths_dtype is float and directed and assortative and init_affinity_filename:
            # case 7: directed + assortative + w from file
            c_v.resize(nof_vertices, nof_groups)
            report.c_obj = c_multitensor_factorization[
                bidirectionalS,
                DiagonalTensor[numpy.float_t],
                init_symmetric_tensor_from_initial[DiagonalTensor[numpy.float_t]],
                vertex_t,
                numpy.float_t
            ](
                < const vector[vertex_t] & > edges_start,
                < const vector[vertex_t] & > edges_end,
                < const vector[numpy.float_t] & > edges_weights,
                nof_realizations, max_nof_iterations, nof_convergences,
                labels, c_u, c_v, c_affinity, deref(rng)
            )
    finally:
        # delete pointers
        del rng

    # U and V outputs
    u = numpy.array(
        [c_u(i, j) for i in range(c_u.get_nrows()) for j in range(c_u.get_ncols())]
    ).reshape((c_u.get_nrows(), c_u.get_ncols()))
    v = numpy.array(
        [c_v(i, j) for i in range(c_v.get_nrows()) for j in range(c_v.get_ncols())]
    ).reshape((c_v.get_nrows(), c_v.get_ncols()))

    # Affinity output
    # We return a format similar to that returned
    # by the function read_affinity_data
    # i.e. a list of arrays
    affinity_ravel = numpy.array(
        c_affinity
    )
    num_vals = affinity_ravel.size // nof_layers
    affinity = []
    for l in range(nof_layers):
        begin = l * num_vals
        end = (l + 1) * num_vals
        affinity.append(
            numpy.array(affinity_ravel[begin:end]).reshape((-1, nof_groups))
        )
    return u, v, affinity, report

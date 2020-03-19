from libcpp.string cimport string
from libcpp.vector cimport vector

import numpy
cimport numpy

from .Report cimport Report
from .tensor cimport Tensor, Matrix, SymmetricTensor, DiagonalTensor


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


# This is a python wrapper for the report returned by the solver. For
# a description of the C interface you can refer to Report.pxd.
cdef class PyReport:
    cdef Report report

    @property
    def nof_realizations(self):
        return self.report.nof_realizations

    @nof_realizations.setter
    def nof_realizations(self, nof_realizations):
        self.report.nof_realizations = nof_realizations

    @property
    def vec_L2(self):
        return self.report.vec_L2

    @vec_L2.setter
    def vec_L2(self, vec_L2):
        self.report.vec_L2 = vec_L2

    @property
    def vec_iter(self):
        return self.report.vec_iter

    @vec_iter.setter
    def vec_iter(self, vec_iter):
        self.report.vec_iter = vec_iter

    @property
    def vec_term_reason(self):
        return self.report.vec_term_reason

    @vec_term_reason.setter
    def vec_term_reason(self, vec_term_reason):
        self.report.vec_term_reason = vec_term_reason

    @property
    def duration(self):
        return self.report.duration

    @duration.setter
    def duration(self, duration):
        self.report.duration = duration

    def max_L2(self):
        return self.report.max_L2()


# A small utility function for calculating the number of network
# vertices from the edge endpoints.
cdef extern from "multitensor/utils.hpp" namespace "multitensor::utils":
    cdef size_t get_num_vertices[vertex_t](
        const vector[vertex_t] &edges_start,
        const vector[vertex_t] &edges_end
    )


# Main entry point of the algorithm.
cdef extern from "multitensor/main.hpp" namespace "multitensor":
    cdef Report multitensor_factorization[direction_t, affinity_t, affinity_init_t, vertex_t, weight_t] (
        const vector[vertex_t] &edges_start,
        const vector[vertex_t] &edges_end,
        const vector[weight_t] &edges_weight,
        const size_t &nof_realizations,
        const size_t &max_nof_iterations,
        const size_t &nof_convergences,
        vector[vertex_t] &labels,
        Matrix[numpy.float_t] &u,
        Matrix[numpy.float_t] &v,
        vector[numpy.float_t] &affinity
    )


# Python wrapper for the algorithm.
def py_multitensor_factorization(
    numpy.ndarray[vertex_t, ndim=1, cast=True] edges_start not None,
    numpy.ndarray[vertex_t, ndim=1, cast=True] edges_end not None,
    numpy.ndarray[weight_t, ndim=1, cast=True] edges_weights not None,
    size_t nof_groups,
    numpy.ndarray[numpy.float_t, ndim=1, cast=True] init_affinity,
    bint directed,
    bint assortative,
    size_t nof_realizations,
    size_t max_nof_iterations,
    size_t nof_convergences,
):

    # We start by initializing the output variables of the C function.
    # Report can be initialized as empty.
    cdef Report report

    # The underlying C function deduces the number of groups from the
    # shane of the affinity matrix passed as an input. For a python
    # version we want to keep a more python approach and pass an
    # explicit parameter. This implies that we have to repeat the math
    # done in C function in reverse.
    cdef size_t nof_edges = edges_start.size
    cdef size_t nof_vertices = get_num_vertices[vertex_t](
            <const vector[vertex_t] &> edges_start,
            <const vector[vertex_t] &> edges_end)
    cdef size_t nof_layers = edges_weights.size / nof_edges

    # Preallocate the output matrices.
    cdef Matrix[numpy.float_t] u = Matrix[numpy.float_t](nof_vertices, nof_groups)
    cdef Matrix[numpy.float_t] v = Matrix[numpy.float_t](0, 0)

    # Preallocate the affinity vector depending on the
    # assortative/nonassortative case.
    cdef vector[numpy.float_t] affinity

    if assortative:
        affinity_size = nof_groups * nof_layers
    else:
        affinity_size = nof_groups * nof_groups * nof_layers

    # Initialize the affinity vector
    if init_affinity is None:
        # as an empty one if there is no initial condition
        affinity = vector[numpy.float_t](affinity_size)
    else:
        # or by copying the initial condition vector
        affinity = <vector[numpy.float_t]> init_affinity

    # Create an empty label vector. Not sure why we need this :)
    cdef vector[vertex_t] labels = vector[vertex_t](nof_vertices)

    # The code below reproduces the swithc case from MultiTensor.hpp
    # (starting on line 175 as of moment of writing). This is done for
    # the case of weight_t being an int.
    if weight_t is numpy.int_t and not directed and not assortative and init_affinity is None:
        # case 0: undirected + non-assortative + w random
        report = multitensor_factorization[
            undirectedS,
            SymmetricTensor[numpy.float_t],
            init_symmetric_tensor_random,
            vertex_t,
            numpy.int_t
        ](
            <const vector[vertex_t] &> edges_start,
            <const vector[vertex_t] &> edges_end,
            <const vector[weight_t] &> edges_weights,
            nof_realizations, max_nof_iterations, nof_convergences,
            labels, u, v, affinity
        )

    if weight_t is numpy.int_t and directed and not assortative and init_affinity is None:
        # case 1: directed + non-assortative + w random
        v.resize(nof_vertices, nof_groups)
        report = multitensor_factorization[
            bidirectionalS,
            SymmetricTensor[numpy.float_t],
            init_symmetric_tensor_random,
            vertex_t,
            numpy.int_t
        ](
            <const vector[vertex_t] &> edges_start,
            <const vector[vertex_t] &> edges_end,
            <const vector[weight_t] &> edges_weights,
            nof_realizations, max_nof_iterations, nof_convergences,
            labels, u, v, affinity
        )

    if weight_t is numpy.int_t and not directed and assortative and init_affinity is None:
        # case 2: undirected + assortative + w random
        report = multitensor_factorization[
            undirectedS,
            DiagonalTensor[numpy.float_t],
            init_symmetric_tensor_random,
            vertex_t,
            numpy.int_t
        ](
            <const vector[vertex_t] &> edges_start,
            <const vector[vertex_t] &> edges_end,
            <const vector[weight_t] &> edges_weights,
            nof_realizations, max_nof_iterations, nof_convergences,
            labels, u, v, affinity
        )

    if weight_t is numpy.int_t and directed and assortative and init_affinity is None:
        # case 3: directed + assortative + w random
        v.resize(nof_vertices, nof_groups)
        report = multitensor_factorization[
            bidirectionalS,
            DiagonalTensor[numpy.float_t],
            init_symmetric_tensor_random,
            vertex_t,
            numpy.int_t
        ](
            <const vector[vertex_t] &> edges_start,
            <const vector[vertex_t] &> edges_end,
            <const vector[weight_t] &> edges_weights,
            nof_realizations, max_nof_iterations, nof_convergences,
            labels, u, v, affinity
        )

    if weight_t is numpy.int_t and not directed and not assortative and init_affinity is not None:
        # case 4: undirected + non-assortative + w from file
        report = multitensor_factorization[
            undirectedS,
            SymmetricTensor[numpy.float_t],
            init_symmetric_tensor_from_initial[SymmetricTensor[numpy.float_t]],
            vertex_t,
            numpy.int_t
        ](
            <const vector[vertex_t] &> edges_start,
            <const vector[vertex_t] &> edges_end,
            <const vector[weight_t] &> edges_weights,
            nof_realizations, max_nof_iterations, nof_convergences,
            labels, u, v, affinity
        )

    if weight_t is numpy.int_t and directed and not assortative and init_affinity is not None:
        # case 5: directed + non-assortative + w from file
        v.resize(nof_vertices, nof_groups)
        report = multitensor_factorization[
            bidirectionalS,
            SymmetricTensor[numpy.float_t],
            init_symmetric_tensor_from_initial[SymmetricTensor[numpy.float_t]],
            vertex_t,
            numpy.int_t
        ](
            <const vector[vertex_t] &> edges_start,
            <const vector[vertex_t] &> edges_end,
            <const vector[weight_t] &> edges_weights,
            nof_realizations, max_nof_iterations, nof_convergences,
            labels, u, v, affinity
        )

    if weight_t is numpy.int_t and not directed and assortative and init_affinity is not None:
        # case 6: undirected + assortative + w from file
        report = multitensor_factorization[
            undirectedS,
            DiagonalTensor[numpy.float_t],
            init_symmetric_tensor_from_initial[DiagonalTensor[numpy.float_t]],
            vertex_t,
            numpy.int_t
        ](
            <const vector[vertex_t] &> edges_start,
            <const vector[vertex_t] &> edges_end,
            <const vector[weight_t] &> edges_weights,
            nof_realizations, max_nof_iterations, nof_convergences,
            labels, u, v, affinity
        )

    if weight_t is numpy.int_t and directed and assortative and init_affinity is not None:
        # case 7: directed + assortative + w from file
        v.resize(nof_vertices, nof_groups)
        report = multitensor_factorization[
            bidirectionalS,
            DiagonalTensor[numpy.float_t],
            init_symmetric_tensor_from_initial[DiagonalTensor[numpy.float_t]],
            vertex_t,
            numpy.int_t
        ](
            <const vector[vertex_t] &> edges_start,
            <const vector[vertex_t] &> edges_end,
            <const vector[weight_t] &> edges_weights,
            nof_realizations, max_nof_iterations, nof_convergences,
            labels, u, v, affinity
        )

    # At the moment Cython does not support using fused types to
    # instantiate a template argument, and we have to explicitly pick
    # the implementation for every type variant. This is the same
    # switch case repeated for floating-piont weight_t.
    if weight_t is numpy.float_t and not directed and not assortative and init_affinity is None:
        # case 0: undirected + non-assortative + w random
        report = multitensor_factorization[
            undirectedS,
            SymmetricTensor[numpy.float_t],
            init_symmetric_tensor_random,
            vertex_t,
            numpy.float_t
        ](
            <const vector[vertex_t] &> edges_start,
            <const vector[vertex_t] &> edges_end,
            <const vector[weight_t] &> edges_weights,
            nof_realizations, max_nof_iterations, nof_convergences,
            labels, u, v, affinity
        )

    if weight_t is numpy.float_t and directed and not assortative and init_affinity is None:
        # case 1: directed + non-assortative + w random
        v.resize(nof_vertices, nof_groups)
        report = multitensor_factorization[
            bidirectionalS,
            SymmetricTensor[numpy.float_t],
            init_symmetric_tensor_random,
            vertex_t,
            numpy.float_t
        ](
            <const vector[vertex_t] &> edges_start,
            <const vector[vertex_t] &> edges_end,
            <const vector[weight_t] &> edges_weights,
            nof_realizations, max_nof_iterations, nof_convergences,
            labels, u, v, affinity
        )

    if weight_t is numpy.float_t and not directed and assortative and init_affinity is None:
        # case 2: undirected + assortative + w random
        report = multitensor_factorization[
            undirectedS,
            DiagonalTensor[numpy.float_t],
            init_symmetric_tensor_random,
            vertex_t,
            numpy.float_t
        ](
            <const vector[vertex_t] &> edges_start,
            <const vector[vertex_t] &> edges_end,
            <const vector[weight_t] &> edges_weights,
            nof_realizations, max_nof_iterations, nof_convergences,
            labels, u, v, affinity
        )

    if weight_t is numpy.float_t and directed and assortative and init_affinity is None:
        # case 3: directed + assortative + w random
        v.resize(nof_vertices, nof_groups)
        report = multitensor_factorization[
            bidirectionalS,
            DiagonalTensor[numpy.float_t],
            init_symmetric_tensor_random,
            vertex_t,
            numpy.float_t
        ](
            <const vector[vertex_t] &> edges_start,
            <const vector[vertex_t] &> edges_end,
            <const vector[weight_t] &> edges_weights,
            nof_realizations, max_nof_iterations, nof_convergences,
            labels, u, v, affinity
        )

    if weight_t is numpy.float_t and not directed and not assortative and init_affinity is not None:
        # case 4: undirected + non-assortative + w from file
        report = multitensor_factorization[
            undirectedS,
            SymmetricTensor[numpy.float_t],
            init_symmetric_tensor_from_initial[SymmetricTensor[numpy.float_t]],
            vertex_t,
            numpy.float_t
        ](
            <const vector[vertex_t] &> edges_start,
            <const vector[vertex_t] &> edges_end,
            <const vector[weight_t] &> edges_weights,
            nof_realizations, max_nof_iterations, nof_convergences,
            labels, u, v, affinity
        )

    if weight_t is numpy.float_t and directed and not assortative and init_affinity is not None:
        # case 5: directed + non-assortative + w from file
        v.resize(nof_vertices, nof_groups)
        report = multitensor_factorization[
            bidirectionalS,
            SymmetricTensor[numpy.float_t],
            init_symmetric_tensor_from_initial[SymmetricTensor[numpy.float_t]],
            vertex_t,
            numpy.float_t
        ](
            <const vector[vertex_t] &> edges_start,
            <const vector[vertex_t] &> edges_end,
            <const vector[weight_t] &> edges_weights,
            nof_realizations, max_nof_iterations, nof_convergences,
            labels, u, v, affinity
        )

    if weight_t is numpy.float_t and not directed and assortative and init_affinity is not None:
        # case 6: undirected + assortative + w from file
        report = multitensor_factorization[
            undirectedS,
            DiagonalTensor[numpy.float_t],
            init_symmetric_tensor_from_initial[DiagonalTensor[numpy.float_t]],
            vertex_t,
            numpy.float_t
        ](
            <const vector[vertex_t] &> edges_start,
            <const vector[vertex_t] &> edges_end,
            <const vector[weight_t] &> edges_weights,
            nof_realizations, max_nof_iterations, nof_convergences,
            labels, u, v, affinity
        )

    if weight_t is numpy.float_t and directed and assortative and init_affinity is not None:
        # case 7: directed + assortative + w from file
        v.resize(nof_vertices, nof_groups)
        report = multitensor_factorization[
            bidirectionalS,
            DiagonalTensor[numpy.float_t],
            init_symmetric_tensor_from_initial[DiagonalTensor[numpy.float_t]],
            vertex_t,
            numpy.float_t
        ](
            <const vector[vertex_t] &> edges_start,
            <const vector[vertex_t] &> edges_end,
            <const vector[weight_t] &> edges_weights,
            nof_realizations, max_nof_iterations, nof_convergences,
            labels, u, v, affinity
        )

    # Finally, we wrap all the return arguments as proper python
    # objects. It is pretty much the same for every data type in here
    # -- we initialize an empty python instance and then we manually
    # copy all the data from C object into a python object.
    py_report = PyReport()
    py_report.nof_realizations = report.nof_realizations
    py_report.vec_L2 = report.vec_L2
    py_report.vec_iter = report.vec_iter
    py_report.vec_term_reason = report.vec_term_reason
    py_report.duration = report.duration

    cdef numpy.ndarray py_u = numpy.zeros(
        [u.get_nrows(), u.get_ncols()], dtype=numpy.float)

    cdef numpy.ndarray py_v = numpy.zeros(
        [v.get_nrows(), v.get_ncols()], dtype=numpy.float)

    cdef numpy.ndarray py_affinity = numpy.zeros(
        [affinity.size()], dtype=numpy.float)

    for i in range(u.get_nrows()):
        for j in range(u.get_ncols()):
            py_u[i, j] = u(i, j)

    for i in range(v.get_nrows()):
        for j in range(v.get_ncols()):
            py_v[i, j] = v(i, j)

    for i in range(affinity.size()):
        py_affinity[i] = affinity[i]

    return py_u, py_v, py_affinity, py_report

import numpy

from .main import py_multitensor_factorization, PyReport as Report
from .app_utils import py_read_adjacency_data, py_read_affinity_data


__all__ = ["multitensor_factorization", "read_adjacency_data", "read_affinity_data"]


def read_adjacency_data(filename):
    return py_read_adjacency_data(filename.encode("utf-8"))


def read_affinity_data(filename, num_groups, num_layers, assortative=False):
    return py_read_affinity_data(
        filename.encode("utf-8"),
        num_groups, num_layers, assortative)


def multitensor_factorization(
    edges_start,
    edges_end,
    edges_weights,
    nof_groups,
    init_affinity=None,
    directed=False,
    assortative=False,
    nof_realizations=1,
    max_nof_iterations=500,
    nof_convergences=10
):

    edges_weights = numpy.ravel(edges_weights)

    return py_multitensor_factorization(
        edges_start,
        edges_end,
        edges_weights,
        nof_groups,
        init_affinity,
        directed,
        assortative,
        nof_realizations,
        max_nof_iterations,
        nof_convergences)

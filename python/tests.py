# Copyright (c) 2019, Max Planck Society / Software Workshop - Max Planck Institute for Intelligent Systems
# Distributed under the GNU GPL license version 3
# See file LICENSE.md or at https://github.com/MPI-IS/multitensor/LICENSE.md

"""Module for testing the Python extension."""

from pathlib import Path
from unittest import TestCase

import numpy

from multitensor import (
    multitensor_factorization,
    read_adjacency_data,
    read_affinity_data,
    ReportWrapper
)

PARENT_DIRECTORY = Path(__file__).parent.parent
DATA_DIRECTORY = PARENT_DIRECTORY / "data"


class InputFileMixin:
    """
    This is a bare-bones test case using a set of given input
    files. Though the canonical output files are given as well, at the
    moment we don't do any comparison with them, since writing data
    files is hard :) At the moment the only thing we check is that the
    algorithm works for those datafiles and that the output has the
    expected type and the shape.
    """

    # Attributes required in the derived classes
    adjacency_filename = None
    num_groups = None

    # Optinal arguments
    affinity_filename = None
    assortative = False
    undirected = False
    num_realizations = 1
    num_iterations = 100
    num_convergences = 10

    def setUp(self):

        adjacency_filename = DATA_DIRECTORY / self.adjacency_filename
        self.edges_start, self.edges_end, self.edges_weights = (
            read_adjacency_data(adjacency_filename))

        self.num_vertices = len(
            set(self.edges_start) | set(self.edges_end))

        self.num_layers = int(len(self.edges_weights) / len(self.edges_start))

        self.init_affinity = None
        if self.affinity_filename:
            affinity_filename = DATA_DIRECTORY / self.affinity_filename
            self.init_affinity = read_affinity_data(
                affinity_filename, self.num_groups, self.num_layers, self.assortative)

    def test_factorization_does_not_crash_and_has_sensible_output(self):
        """
        Run the algorithm on a given input, wait for it to finish without
        crashing and then verify that everything has a reasonable type
        and shape.
        """
        u, v, w, report = multitensor_factorization(
            self.edges_start,
            self.edges_end,
            self.edges_weights,
            self.num_groups,
            not self.undirected,
            self.assortative,
            self.num_realizations,
            self.num_iterations,
            self.num_convergences,
            self.init_affinity)

        self.assertIsInstance(u, numpy.ndarray)
        self.assertIsInstance(v, numpy.ndarray)
        self.assertIsInstance(w, numpy.ndarray)
        self.assertIsInstance(report, ReportWrapper)

        self.assertEqual(u.shape, (self.num_vertices, self.num_groups))
        if not self.undirected:
            self.assertEqual(v.shape, (self.num_vertices, self.num_groups))
        if self.assortative:
            self.assertEqual(w.shape, (self.num_groups * self.num_layers, ))
        else:
            self.assertEqual(w.shape, (self.num_groups * self.num_groups * self.num_layers, ))


class DirectedMainInputFileTestCase(InputFileMixin, TestCase):
    adjacency_filename = "main/adjacency.dat"
    num_groups = 2


class UndirectedInputFileTestCase(InputFileMixin, TestCase):
    adjacency_filename = "undirected/adjacency.dat"
    num_groups = 2
    undirected = True


class WInputInputFileTestCase(InputFileMixin, TestCase):
    adjacency_filename = "w_input/adjacency_k2L4.dat"
    affinity_filename = "w_input/w_k2_k2L4_A.dat"
    num_groups = 2


class MultiRealInputFileTestCase(InputFileMixin, TestCase):
    adjacency_filename = "multi_real/adjacency_k2L4.dat"
    affinity_filename = "multi_real/w_k2_k2L4_r2.dat"
    num_groups = 2
    num_realizations = 2


class AssortativeInputFileTestCase(InputFileMixin, TestCase):
    adjacency_filename = "assortative/adjacency_assortative_k3L4.dat"
    num_groups = 3
    num_realizations = 1
    assortative = True

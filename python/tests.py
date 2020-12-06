# Copyright (c) 2019, Max Planck Society / Software Workshop - Max Planck Institute for Intelligent Systems
# Distributed under the GNU GPL license version 3
# See file LICENSE.md or at https://github.com/MPI-IS/multitensor/LICENSE.md
"""Module for testing the Python extension."""

import os
from pathlib import Path
from unittest import TestCase

import numpy as np

from multitensor import ReportWrapper, run

DATA_DIRECTORY = Path(os.environ['ROOT_DIR']) / "data"
TEST_DATA_DIRECTORY = Path(os.environ['ROOT_DIR']) / "tests" / "data"
RELATIVE_TOL = 1e-5


class InputFileMixin:
    """
    This is a bare - bones test case using a set of given input
    files. Though the canonical output files are given as well, at the
    moment we don't do any comparison with them, since writing data
    files is hard: ) At the moment the only thing we check is that the
    algorithm works for those datafiles and that the output has the
    expected type and the shape.
    """

    # Attributes required in the derived classes
    adjacency_filename = None
    num_groups = None

    # Optinal arguments
    affinity_filename = None
    assortative = False
    directed = True
    num_realizations = 1
    num_iterations = 100
    seed = 5489

    # Data for comparison with C++ runs
    u_compare_filename = None
    v_compare_filename = None
    w_compare_filename = None

    def test_factorization_does_not_crash_and_has_sensible_output(self):
        """
        Run the algorithm on a given input, wait for it to finish without
        crashing and then verify that everything has a reasonable type
        and shape.
        """
        # Initial affinity data
        w_file = None
        if self.affinity_filename:
            w_file = DATA_DIRECTORY / self.folder / self.affinity_filename

        u, v, w, report = run(
            DATA_DIRECTORY / self.folder / self.adjacency_filename,
            self.num_groups,
            self.directed,
            self.assortative,
            self.num_realizations,
            self.num_iterations,
            init_affinity_filename=w_file,
            seed=self.seed)

        # Check report
        self.assertIsInstance(report, ReportWrapper)
        self.assertEqual(report.nof_realizations, self.num_realizations)
        self.assertEqual(len(report.vec_iter), report.nof_realizations)
        for it in report.vec_iter:
            self.assertLessEqual(it, self.num_iterations)

        # Check u and w
        self.assertIsInstance(u, np.ndarray)
        num_vertices = u.shape[0]
        num_groups = u.shape[1] - 1
        self.assertEqual(num_groups, self.num_groups)

        self.assertIsInstance(w, list)
        self.assertIsInstance(w[0], np.ndarray)
        for l in w:
            if self.assortative:
                self.assertEqual(l.shape, (num_groups,))
            else:
                self.assertEqual(l.shape, (num_groups, num_groups))

        u_compare = np.loadtxt(TEST_DATA_DIRECTORY / self.folder / self.u_compare_filename)
        w_compare = np.loadtxt(
            TEST_DATA_DIRECTORY / self.folder / self.w_compare_filename, comments=['#', "a="]
        )
        self.assertTrue((u[:, 0] == u_compare[:, 0]).all())
        np.testing.assert_allclose(u[:, 1:], u_compare[:, 1:], rtol=RELATIVE_TOL)
        for l, w_layer in enumerate(w):
            np.testing.assert_allclose(
                w_layer, w_compare[num_groups * l:num_groups * (l + 1)], rtol=RELATIVE_TOL)

        # Check v
        if self.directed:
            self.assertIsInstance(v, np.ndarray)
            self.assertEqual(v.shape, (num_vertices, num_groups + 1))
            v_compare = np.loadtxt(TEST_DATA_DIRECTORY / self.folder / self.v_compare_filename)
            self.assertTrue((v[:, 0] == v_compare[:, 0]).all())
            np.testing.assert_allclose(v[:, 1:], v_compare[:, 1:], rtol=RELATIVE_TOL)
        else:
            self.assertIsNone(v)


class TestDirectedMainInputFile(InputFileMixin, TestCase):
    folder = "main"
    adjacency_filename = "adjacency.dat"
    num_groups = 2
    u_compare_filename = "u_K2_compare.txt"
    v_compare_filename = "v_K2_compare.txt"
    w_compare_filename = "w_K2_compare.txt"


class TestUndirectedInputFile(InputFileMixin, TestCase):
    folder = "undirected"
    adjacency_filename = "adjacency.dat"
    num_groups = 2
    directed = False
    u_compare_filename = "u_K2_compare.txt"
    w_compare_filename = "w_K2_compare.txt"


class TestWInputInputFile(InputFileMixin, TestCase):
    folder = "w_input"
    adjacency_filename = "adjacency_k2L4.dat"
    affinity_filename = "w_k2_k2L4_A.dat"
    num_groups = 2
    u_compare_filename = "u_K2_k2L4_i2_compare.txt"
    v_compare_filename = "v_K2_k2L4_i2_compare.txt"
    w_compare_filename = "w_K2_k2L4_i2_compare.txt"


class TestMultiRealInputFile(InputFileMixin, TestCase):
    folder = "multi_real"
    adjacency_filename = "adjacency_k2L4.dat"
    affinity_filename = "w_k2_k2L4_r2.dat"
    num_groups = 2
    num_realizations = 2
    u_compare_filename = "u_K2_k2L4_r2_i2_compare.txt"
    v_compare_filename = "v_K2_k2L4_r2_i2_compare.txt"
    w_compare_filename = "w_K2_k2L4_r2_i2_compare.txt"


class TestAssortativeInputFile(InputFileMixin, TestCase):
    folder = "assortative"
    adjacency_filename = "adjacency_assortative_k3L4.dat"
    num_groups = 3
    num_realizations = 1
    assortative = True
    u_compare_filename = "u_K3_adjass_k3L4_compare.txt"
    v_compare_filename = "v_K3_adjass_k3L4_compare.txt"
    w_compare_filename = "w_K3_adjass_k3L4_compare.txt"


class TestErrorHandling(TestCase):
    """Class for testing how errors are handled."""

    folder = "main"
    adjacency_filename = "adjacency.dat"
    num_groups = 1  # that should cause an error

    def test_error(self):
        """Checks that the appropriate exception is raised."""

        with self.assertRaises(RuntimeError):
            _ = run(DATA_DIRECTORY / self.folder / self.adjacency_filename, self.num_groups)

# Copyright (c) 2019, Max Planck Society / Software Workshop - Max Planck Institute for Intelligent Systems
# Distributed under the GNU GPL license version 3
# See file LICENSE.md or at https://github.com/MPI-IS/multitensor/LICENSE.md

"""
Module for running and analyzing a benchmark.

:author: Jean-Claude Passy <jean-claude.passy@tuebingen.mpg.de>
"""


import argparse
from itertools import combinations
import logging
from pathlib import Path
import subprocess

import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format='*** [%(levelname)s] %(asctime)s %(message)s'
)


def cosine_similarity(U_infer, U0, P):
    """Function calculation the cosine similarity."""

    # Dimensions
    N, K = U0.shape

    # Permute infered matrix
    # TODO: should use P instead of recalculating a permutation
    P2 = calculate_permuation(U_infer, U0)
    U_infer = U_infer @ P2

    # Normalize using the infinite norm
    norm_inf = np.linalg.norm(U_infer, axis=1)
    norm0 = np.linalg.norm(U0, axis=1)
    U_infer = U_infer / norm_inf.reshape(U_infer.shape[0], 1)
    U0 = U0 / norm0.reshape(U0.shape[0], 1)

    # Cosine similarity
    CS = sum(
        (U_infer[:, k].T @ U0[:, k] for k in range(K))
    )

    # OLD version
    # CS1 = sum(
    #     (np.dot(U_infer[i],U0[i]) for i in range(N) )
    # )
    # assert abs(CS/N-CS1/N)<0.01
    return CS / N


def calculate_permuation(U_infer, U0):
    """Function calculating the permutation to match the groups."""

    # Dimensions
    N, RANK = U0.shape
    M = U_infer.T @ U0 / float(N)  # dim=RANKxRANK
    rows = np.zeros(RANK)
    columns = np.zeros(RANK)
    P = np.zeros((RANK, RANK))  # Permutation matrix
    for t in range(RANK):
        # Find the max element in the remaining submatrix,
        # the one with rows and columns removed from previous iterations
        max_entry = 0.
        c_index = 1
        r_index = 1
        for i in range(RANK):
            if columns[i] == 0:
                for j in range(RANK):
                    if rows[j] == 0:
                        if M[j, i] > max_entry:
                            max_entry = M[j, i]
                            c_index = i
                            r_index = j

        P[r_index, c_index] = 1
        columns[c_index] = 1
        rows[r_index] = 1

    return P


class BenchmarkTest:
    """Class for a benchmark run."""

    UNAME = "u"
    VNAME = "v"
    WNAME = "w"
    COMPONENTS = (UNAME, VNAME, WNAME)

    ADJACENCY_FILENAME_PREFIX = 'adjacency'
    PATTERN_FILENAME = '_K{k}_{run}_dis_nn.dat'
    RESULTS_FOLDER = "results"

    def __init__(self, exe, num_groups, num_layers,
                 num_nodes=300, num_real=10, seed=5489,
                 name='BenchmarkTest', directory=None):

        self.exe = Path(exe)
        self.num_groups = num_groups
        self.num_layers = num_layers
        self.num_nodes = num_nodes
        self.num_real = num_real
        self.seed = seed
        self.name = name
        self.directory = Path(directory) if directory else Path.cwd()

        # Check executable exists
        if not self.exe.is_file():
            raise RuntimeError(f"Cannot find executable {self.exe}")

        # Calculate sample size - number of adjacency files
        self.adjacency_files = sorted(
            self.directory.glob(self.ADJACENCY_FILENAME_PREFIX + "*")
        )
        self.sample_size = len(self.adjacency_files)
        if not self.sample_size:
            raise RuntimeError(
                "Cannot find adjacency files "
                f"with pattern {self.ADJACENCY_FILENAME_PREFIX} "
                f"in {self.directory}"
            )

        # Build output directory if necessary
        self.results_path = self.directory / self.RESULTS_FOLDER
        self.results_path.mkdir(parents=True, exist_ok=True)

    def analyze_results(self):
        nodes = tuple(range(self.num_nodes))
        CS_list = []
        L1_list = []

        for run1, run2 in combinations(range(self.sample_size), 2):

            # Read U memberships
            Us = []
            for run in (run1, run2):
                ufile = self.results_path / (
                    self.UNAME + self.PATTERN_FILENAME.format(k=self.num_groups, run=run)
                )
                Us.append(np.loadtxt(ufile))

            # Re-order by node and remove a dimension
            Us = [u[u[..., 0].argsort()][..., 1:] for u in Us]

            # Permute
            P = calculate_permuation(*Us)
            Us = [Us[0] @ P, Us[1]]

            # Normalize
            Us = [u / u.sum(axis=1).reshape(u.shape[0], 1) for u in Us]

            # Metrics
            CS = cosine_similarity(*(Us + [P]))
            L1 = np.abs(Us[0] - Us[1]).sum() / Us[0].size

            # Save results
            CS_list.append(CS)
            L1_list.append(L1)

        return CS_list, L1_list

    def run(self):
        """Run benchmark test."""

        # Run algorithm on each adjacency file
        # Example of command:
        # Multitensor --k 2 --a adjacency.dat --r 10 --o results
        for i, adj in enumerate(self.adjacency_files):
            logging.info(f"... running with {adj.name}")
            command = [
                self.exe,
                "--k", str(self.num_groups),
                "--a", adj.resolve(),
                "--r", str(self.num_real),
                "--o", self.results_path.resolve(),
                "--s", str(self.seed),
            ]

            subprocess.check_call(command, stdout=subprocess.DEVNULL)

            # Rename output file for simplicity
            for component in self.COMPONENTS:
                new_name = component + self.PATTERN_FILENAME.format(k=self.num_groups, run=i)
                list(self.results_path.glob(component + "_out.dat"))[0].rename(
                    self.results_path / new_name)

        # Analyze results
        CS, L1 = self.analyze_results()

        # Print final results
        S = np.mean(CS)
        sigma_s = np.std(CS)
        L = np.mean(L1)
        sigma_L = np.std(L1)

        results = (
            f"- Results: \n"
            f"\t CS = {S:.4f} +/- {sigma_s:.5f}\n"
            f"\t L1 = {L:.4f} +/- {sigma_L:.5f}\n"
        )
        logging.info(results)


class BenchmarkTestG1(BenchmarkTest):
    """Class for the benchmark with network G1."""

    def __init__(self, exe):
        name = "G1"
        directory = "K2L2"
        num_groups = 2
        num_layers = 2
        super().__init__(exe, num_groups, num_layers, name=name, directory=directory)


class BenchmarkTestG2(BenchmarkTest):
    """Class for the benchmark with network G2."""

    def __init__(self, exe):
        name = "G2"
        directory = "K2L4"
        num_groups = 2
        num_layers = 4
        super().__init__(exe, num_groups, num_layers, name=name, directory=directory)


class BenchmarkTestG3(BenchmarkTest):
    """Class for the benchmark with network G3."""

    def __init__(self, exe):
        name = "G3"
        directory = "K2L4two"
        num_groups = 2
        num_layers = 4
        super().__init__(exe, num_groups, num_layers, name=name, directory=directory)


# Different types of benchmark
all_benchmark_types = {
    'G1': BenchmarkTestG1,
    'G2': BenchmarkTestG2,
    'G3': BenchmarkTestG3,
}


def build_benchmark(*args, **kwargs):
    """Build a benchmark to run."""

    # List of tests to run
    benchmark = []

    # Type of benchmark requested
    test_list = kwargs.pop('test')

    if 'all' in test_list:
        logging.info("All tests will be run.")
        test_list = all_benchmark_types
    for test in test_list:
        benchmark.append(all_benchmark_types[test](**kwargs))

    return benchmark


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Running benchmark",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--exe', metavar='EXE', type=str,
                        default=str(Path.cwd().parent / 'Multitensor'), help='Executable')
    parser.add_argument('--test', action='append', help='Benchmark to run',
                        default=['all'], choices=list(all_benchmark_types.keys()) + ['all'])

    # Get args from command line
    args = parser.parse_args()

    # Build instances to run
    benchmark = build_benchmark(**vars(args))

    # Run the tests
    for test_to_run in benchmark:
        logging.info("=====================================")
        logging.info(f"- Running test {test_to_run.name} from folder {test_to_run.directory}")
        test_to_run.run()

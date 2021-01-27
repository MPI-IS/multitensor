// Copyright (c) 2019, Max Planck Society / Software Workshop - Max Planck Institute for Intelligent Systems
// Distributed under the GNU GPL license version 3
// See file LICENSE.md or at https://github.com/MPI-IS/multitensor/LICENSE.md

/*!
 * @file
 *
 * @brief Testing the Solver class.
 *
 * @author Jean-Claude Passy (jean-claude.passy@tuebingen.mpg.de)
 */

#include <iostream>
#include <random>
#include <ctime>
#include <vector>
#include <cstddef>
#include <boost/test/unit_test.hpp>

#include "multitensor/params.hpp"
#include "multitensor/initialization.hpp"
#include "multitensor/solver.hpp"
#include "fixtures.hpp"

// Fixture
struct fixture_solver : fixture_global
{
    size_t nof_realizations, max_nof_iterations, nof_convergences, nof_groups;

    fixture_solver()
        : nof_realizations(rng() % 10 + 1),
          max_nof_iterations(rng() % 10 + 1),
          nof_convergences(rng() % 10 + 1),
          nof_groups(rng() % 10 + 1)
    {
    }
};

using namespace multitensor::graph;
using namespace multitensor::initialization;
using namespace multitensor::solver;
using namespace multitensor::tensor;

BOOST_FIXTURE_TEST_SUITE(tests_solver, fixture_solver)

// Checks the solver initialization
BOOST_AUTO_TEST_CASE(test_solver_init)
{
    Solver solver(nof_realizations, max_nof_iterations, nof_convergences);
    BOOST_TEST(solver.num_real() == nof_realizations);
    BOOST_TEST(solver.max_iter() == max_nof_iterations);
    BOOST_TEST(solver.num_conv() == nof_convergences);
}

// Checks running the solver on an empty network (corner case)
BOOST_AUTO_TEST_CASE(test_solver_run_empty_network)
{
    Solver solver(nof_realizations, max_nof_iterations, nof_convergences);
    std::vector<size_t> vec_empty{};
    SymmetricTensor<double> w(nof_groups, nof_layers);
    Matrix<double> u(nof_vertices, nof_groups), v(nof_vertices, nof_groups);
    Network A(vec_empty, vec_empty, vec_empty);

    BOOST_CHECK_NO_THROW(solver.run<init_symmetric_tensor_random>(
        *u_list, *v_list, A, u, v, w, rng));

    BOOST_CHECK_NO_THROW(solver.run<init_symmetric_tensor_from_initial<SymmetricTensor<double>>>(
        *u_list, *v_list, A, u, v, w, rng));
}

BOOST_AUTO_TEST_SUITE_END()

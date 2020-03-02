/*!
 * @file
 *
 * @brief Testing the Solver class.
 *
 * @author Jean-Claude Passy (jean-claude.passy@tuebingen.mpg.de)
 */

#include <boost/test/unit_test.hpp>
#include <iostream>
#include <random>
#include <ctime>
#include <vector>

#include "multitensor/parameters.hpp"
#include "multitensor/solver.hpp"

using namespace multitensor::solver;

// Fixture
struct fixture_solver
{
    std::mt19937_64 rng;
    std::time_t seed;
    size_t nof_nodes, nof_groups, nof_layers;
    unsigned int nof_realizations, max_nof_iterations, nof_convergences;

    fixture_solver()
        : seed(std::time(nullptr))
    {
        BOOST_TEST_MESSAGE("In fixture, the seed is " << seed);
        rng.seed(seed);
        nof_realizations = static_cast<unsigned int>((rng() % 10)) + 1;
        max_nof_iterations = static_cast<unsigned int>(rng() % 10) + 1;
        nof_convergences = static_cast<unsigned int>(rng() % 10) + 1;
        nof_nodes = size_t(rng() % 10) + 2;
        nof_groups = size_t(rng() % 10) + 2;
        nof_layers = size_t(rng() % 10) + 1;
    }
};

BOOST_FIXTURE_TEST_SUITE(tests_solver, fixture_solver)

// Checks the solver initialization
BOOST_AUTO_TEST_CASE(test_solver_init)
{
    Solver solver(nof_realizations, max_nof_iterations, nof_convergences);
    BOOST_TEST(solver.num_real() == nof_realizations);
    BOOST_TEST(solver.max_iterations() == max_nof_iterations);
}
// Checks running the solver
BOOST_AUTO_TEST_CASE(test_solver_run)
{
    Solver solver(nof_realizations, max_nof_iterations, nof_convergences);
    //BOOST_CHECK_NO_THROW(solver.run(nof_nodes, nof_groups, nof_layers, std::vector<size_t>(0), std::vector<size_t>(0), rng));
}

BOOST_AUTO_TEST_SUITE_END()

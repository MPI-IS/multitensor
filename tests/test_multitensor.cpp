/*!
 * @file
 *
 * @brief Testing the multitensor algorithm.
 *
 * @author Jean-Claude Passy (jean-claude.passy@tuebingen.mpg.de)
 */

#include <iostream>
#include <random>
#include <ctime>
#include <vector>
#include <string>
#include <cstdlib>
#include <boost/test/unit_test.hpp>

#include "multitensor/parameters.hpp"
#include "multitensor/main.hpp"
#include "fixtures.hpp"

// Fixture
struct fixture_multitensor : fixture_global
{
    std::string output_dir;
    size_t nof_groups, nof_realizations, max_nof_iterations, nof_convergences;

    fixture_multitensor()
        : output_dir("results"),
          nof_groups(rng() % 10 + 2),
          nof_realizations(rng() % 10 + 1),
          max_nof_iterations(rng() % 10 + 1),
          nof_convergences(rng() % 10 + 1)
    {
    }
};

using namespace multitensor;

BOOST_FIXTURE_TEST_SUITE(tests_multitensor, fixture_multitensor)

// Checks that bad input data throws errors
BOOST_AUTO_TEST_CASE(test_errors)
{
    tensor::Tensor<double> w(nof_groups, nof_groups, nof_layers);
    tensor::Tensor<double> u(nof_vertices, nof_groups), v(nof_vertices, nof_groups);
    std::vector<size_t> labels;

    // Not enough edges
    BOOST_CHECK_THROW(
        multitensor_factorization(std::vector<unsigned int>{}, edges_end, edges_weight, false,
                                  nof_realizations, max_nof_iterations, nof_convergences,
                                  labels, u, v, w),
        std::runtime_error);

    // Inconsistent edges (less)
    std::vector<unsigned int> edges_end_small(edges_end.begin(), edges_end.end() - 1);
    BOOST_CHECK_THROW(
        multitensor_factorization(edges_start, edges_end_small, edges_weight, false,
                                  nof_realizations, max_nof_iterations, nof_convergences,
                                  labels, u, v, w),
        std::runtime_error);
    // Inconsistent edges (more)
    edges_end_small.push_back(0);
    edges_end_small.push_back(0);
    BOOST_CHECK_THROW(
        multitensor_factorization(edges_start, edges_end_small, edges_weight, false,
                                  nof_realizations, max_nof_iterations, nof_convergences,
                                  labels, u, v, w),
        std::runtime_error);

    // Not enough vertices
    std::vector<unsigned int> one_vertex{0};
    BOOST_TEST(one_vertex.size() == 1);
    BOOST_CHECK_THROW(
        multitensor_factorization(one_vertex, one_vertex, one_vertex, false,
                                  nof_realizations, max_nof_iterations, nof_convergences,
                                  labels, u, v, w),
        std::runtime_error);

    // Not enough layers
    BOOST_CHECK_THROW(
        multitensor_factorization(edges_start, edges_end, std::vector<unsigned int>{}, false,
                                  nof_realizations, max_nof_iterations, nof_convergences,
                                  labels, u, v, w),
        std::runtime_error);

    // Not enough realizations
    BOOST_CHECK_THROW(
        multitensor_factorization(edges_start, edges_end, edges_weight, false,
                                  0, max_nof_iterations, nof_convergences,
                                  labels, u, v, w),
        std::runtime_error);

    // Not enough iterations
    BOOST_CHECK_THROW(
        multitensor_factorization(edges_start, edges_end, edges_weight, false,
                                  nof_realizations, 0, nof_convergences,
                                  labels, u, v, w),
        std::runtime_error);

    // Not enough convergences needed
    BOOST_CHECK_THROW(
        multitensor_factorization(edges_start, edges_end, edges_weight, false,
                                  nof_realizations, max_nof_iterations, 0,
                                  labels, u, v, w),
        std::runtime_error);

    // OK
    BOOST_CHECK_NO_THROW(
        multitensor_factorization(edges_start, edges_end, edges_weight, false,
                                  nof_realizations, max_nof_iterations, nof_convergences,
                                  labels, u, v, w));
}

BOOST_AUTO_TEST_SUITE_END()

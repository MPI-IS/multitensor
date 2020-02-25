/*!
 * @file
 *
 * @brief Testing the multitensor algorithm.
 *
 * @author Jean-Claude Passy (jean-claude.passy@tuebingen.mpg.de)
 */

#include <boost/test/unit_test.hpp>
#include <iostream>
#include <random>
#include <ctime>
#include <vector>

#include "multitensor_parameters.hpp"
#include "multitensor.hpp"

// Fixture
struct fixture_multitensor
{
    std::mt19937_64 rng;
    std::time_t seed;

    std::vector<unsigned int> edges_in, edges_out, edges_weight;
    size_t nof_nodes, nof_layers, nof_edges, nof_groups;
    unsigned int nof_realizations, max_nof_iterations, nof_convergences;

    fixture_multitensor()
        : seed(std::time(nullptr))
    {
        BOOST_TEST_MESSAGE("In fixture, the seed is " << seed);
        rng.seed(seed);

        nof_nodes = size_t(rng() % 10) + 2;
        nof_layers = size_t(rng() % 10) + 1;
        nof_edges = size_t(rng() % 10) + 1;
        nof_groups = size_t(rng() % 10) + 2;
        nof_realizations = static_cast<unsigned int>(rng() % 10) + 1;
        max_nof_iterations = static_cast<unsigned int>(rng() % 10) + 1;
        nof_convergences = static_cast<unsigned int>(rng() % 10) + 1;

        edges_in.reserve(nof_edges);
        edges_out.reserve(nof_edges);
        edges_weight.reserve(nof_edges * nof_layers);
        for (size_t i = 0; i < nof_edges; i++)
        {
            edges_in.emplace_back(static_cast<unsigned int>(rng() % nof_nodes));
            edges_out.emplace_back(static_cast<unsigned int>(rng() % nof_nodes));
            for (size_t alpha = 0; alpha < nof_layers; alpha++)
            {
                edges_weight.emplace_back(static_cast<unsigned int>(rng() % 100));
            }
        }
    }
};

using namespace multitensor;

BOOST_FIXTURE_TEST_SUITE(tests_multitensor, fixture_multitensor)

// Checks that bad input data throws errors
BOOST_AUTO_TEST_CASE(test_errors)
{
    // Not enough nodes
    BOOST_CHECK_THROW(
        multitensor_algo(edges_in, edges_out, edges_weight,
                         1, nof_layers, nof_groups, nof_realizations, max_nof_iterations, nof_convergences),
        std::runtime_error);

    // Not enough layers
    BOOST_CHECK_THROW(
        multitensor_algo(edges_in, edges_out, edges_weight,
                         nof_nodes, 0, nof_groups, nof_realizations, max_nof_iterations, nof_convergences),
        std::runtime_error);

    // Inconsistent edges
    // Less edges
    std::vector<unsigned int> edges_out_small(edges_out.begin(), edges_out.end() - 1);
    BOOST_CHECK_THROW(
        multitensor_algo(edges_in, edges_out_small, edges_weight,
                         nof_nodes, nof_layers, nof_groups, nof_realizations, max_nof_iterations, nof_convergences),
        std::runtime_error);
    // More edges
    edges_out_small.push_back(0);
    edges_out_small.push_back(0);
    BOOST_CHECK_THROW(
        multitensor_algo(edges_in, edges_out_small, edges_weight,
                         nof_nodes, nof_layers, nof_groups, nof_realizations, max_nof_iterations, nof_convergences),
        std::runtime_error);

    // Inconsistent weights
    // Less weights
    std::vector<unsigned int> weights_small(edges_weight.begin(), edges_weight.end() - 1);
    BOOST_CHECK_THROW(
        multitensor_algo(edges_in, edges_out, weights_small,
                         nof_nodes, nof_layers, nof_groups, nof_realizations, max_nof_iterations, nof_convergences),
        std::runtime_error);
    // More weights
    weights_small.push_back(0);
    weights_small.push_back(0);
    BOOST_CHECK_THROW(
        multitensor_algo(edges_in, edges_out, weights_small,
                         nof_nodes, nof_layers, nof_groups, nof_realizations, max_nof_iterations, nof_convergences),
        std::runtime_error);

    // Not enough groups
    BOOST_CHECK_THROW(
        multitensor_algo(edges_in, edges_out, edges_weight,
                         nof_nodes, nof_layers, 1, nof_realizations, max_nof_iterations, nof_convergences),
        std::runtime_error);

    // Not enough realizations
    BOOST_CHECK_THROW(
        multitensor_algo(edges_in, edges_out, edges_weight,
                         nof_nodes, nof_layers, nof_groups, 0, max_nof_iterations, nof_convergences),
        std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()

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

using namespace multitensor;

// Fixture
struct fixture_multitensor : fixture_global
{
    size_t nof_groups, nof_realizations, max_nof_iterations, nof_convergences;
    std::vector<double> w, w_assort;
    tensor::Matrix<double> u, v;
    std::vector<size_t> labels;

    fixture_multitensor()
        : nof_groups(rng() % 3 + 2),
          nof_realizations(1),
          max_nof_iterations(rng() % 10 + 1),
          nof_convergences(rng() % 10 + 1),
          w(nof_groups * nof_groups * nof_layers),
          w_assort(nof_groups * nof_layers),
          u(nof_vertices, nof_groups),
          v(nof_vertices, nof_groups)
    {
    }
};

BOOST_FIXTURE_TEST_SUITE(tests_multitensor, fixture_multitensor)

// Checks that bad input data throws errors
BOOST_AUTO_TEST_CASE(test_errors)
{
    // Not enough edges
    BOOST_CHECK_THROW(
        multitensor_factorization(std::vector<size_t>{}, edges_end, edges_weight,
                                  nof_realizations, max_nof_iterations, nof_convergences,
                                  labels, u, v, w),
        std::runtime_error);

    // Inconsistent edges (less)
    std::vector<size_t> edges_end_small(edges_end.begin(), edges_end.end() - 1);
    BOOST_CHECK_THROW(
        multitensor_factorization(edges_start, edges_end_small, edges_weight,
                                  nof_realizations, max_nof_iterations, nof_convergences,
                                  labels, u, v, w),
        std::runtime_error);
    // Inconsistent edges (more)
    edges_end_small.push_back(0);
    edges_end_small.push_back(0);
    BOOST_CHECK_THROW(
        multitensor_factorization(edges_start, edges_end_small, edges_weight,
                                  nof_realizations, max_nof_iterations, nof_convergences,
                                  labels, u, v, w),
        std::runtime_error);

    // Not enough vertices
    std::vector<size_t> one_vertex{0};
    BOOST_TEST(one_vertex.size() == 1);
    BOOST_CHECK_THROW(
        multitensor_factorization(one_vertex, one_vertex, one_vertex,
                                  nof_realizations, max_nof_iterations, nof_convergences,
                                  labels, u, v, w),
        std::runtime_error);

    // Not enough layers
    BOOST_CHECK_THROW(
        multitensor_factorization(edges_start, edges_end, std::vector<size_t>{},
                                  nof_realizations, max_nof_iterations, nof_convergences,
                                  labels, u, v, w),
        std::runtime_error);

    // Not enough realizations
    BOOST_CHECK_THROW(
        multitensor_factorization(edges_start, edges_end, edges_weight,
                                  0, max_nof_iterations, nof_convergences,
                                  labels, u, v, w),
        std::runtime_error);

    // Not enough iterations
    BOOST_CHECK_THROW(
        multitensor_factorization(edges_start, edges_end, edges_weight,
                                  nof_realizations, 0, nof_convergences,
                                  labels, u, v, w),
        std::runtime_error);

    // Not enough convergences needed
    BOOST_CHECK_THROW(
        multitensor_factorization(edges_start, edges_end, edges_weight,
                                  nof_realizations, max_nof_iterations, 0,
                                  labels, u, v, w),
        std::runtime_error);

    // OK
    BOOST_CHECK_NO_THROW(
        multitensor_factorization(edges_start, edges_end, edges_weight,
                                  nof_realizations, max_nof_iterations, nof_convergences,
                                  labels, u, v, w));
}

// Checks the different algorithm types
BOOST_AUTO_TEST_CASE(test_algo_types)
{
    using namespace boost;
    using namespace multitensor::initialization;
    using namespace multitensor::tensor;

    // To make things faster
    nof_realizations = 2;

    // Undirected + non-assortative + w random
    BOOST_CHECK_NO_THROW(
        multitensor_factorization<undirectedS>(
            edges_start, edges_end, edges_weight,
            nof_realizations, max_nof_iterations, nof_convergences,
            labels, u, v, w));

    // Directed + non-assortative + w random
    BOOST_CHECK_NO_THROW(
        multitensor_factorization(
            edges_start, edges_end, edges_weight,
            nof_realizations, max_nof_iterations, nof_convergences,
            labels, u, v, w));

    // Undirected + assortative + w random
    BOOST_CHECK_NO_THROW(
        (multitensor_factorization<undirectedS,
                                   DiagonalTensor<double>>)(

            edges_start, edges_end, edges_weight,
            nof_realizations, max_nof_iterations, nof_convergences,
            labels, u, v, w_assort));

    // Directed + assortative + w random
    BOOST_CHECK_NO_THROW(
        (multitensor_factorization<bidirectionalS,
                                   DiagonalTensor<double>>)(

            edges_start, edges_end, edges_weight,
            nof_realizations, max_nof_iterations, nof_convergences,
            labels, u, v, w_assort));

    // Undirected + non-assortative + w from file
    BOOST_CHECK_NO_THROW(
        (multitensor_factorization<undirectedS,
                                   SymmetricTensor<double>,
                                   init_symmetric_tensor_from_initial<SymmetricTensor<double>>>)(

            edges_start, edges_end, edges_weight,
            nof_realizations, max_nof_iterations, nof_convergences,
            labels, u, v, w));

    // Directed + non-assortative + w from file
    BOOST_CHECK_NO_THROW(
        (multitensor_factorization<bidirectionalS,
                                   SymmetricTensor<double>,
                                   init_symmetric_tensor_from_initial<SymmetricTensor<double>>>)(

            edges_start, edges_end, edges_weight,
            nof_realizations, max_nof_iterations, nof_convergences,
            labels, u, v, w));

    // Undirected + assortative + w from file
    BOOST_CHECK_NO_THROW(
        (multitensor_factorization<undirectedS,
                                   DiagonalTensor<double>,
                                   init_symmetric_tensor_from_initial<DiagonalTensor<double>>>)(

            edges_start, edges_end, edges_weight,
            nof_realizations, max_nof_iterations, nof_convergences,
            labels, u, v, w_assort));

    // Directed + assortative + w from file
    BOOST_CHECK_NO_THROW(
        (multitensor_factorization<bidirectionalS,
                                   DiagonalTensor<double>,
                                   init_symmetric_tensor_from_initial<DiagonalTensor<double>>>)(

            edges_start, edges_end, edges_weight,
            nof_realizations, max_nof_iterations, nof_convergences,
            labels, u, v, w_assort));
}

BOOST_AUTO_TEST_SUITE_END()

/*!
 * @file
 *
 * @brief Testing the Tensory class.
 *
 * @author Jean-Claude Passy (jean-claude.passy@tuebingen.mpg.de)
 */

#define BOOST_TEST_MODULE tests_multitensor
#include <boost/test/unit_test.hpp>
#include <random>
#include <ctime>

#include "multitensor_parameters.hpp"
#include "multitensor_tensor.hpp"

// Fixture
struct fixture
{
    std::mt19937_64 rng;
    std::time_t seed;

    fixture()
        : seed(std::time(nullptr))
    {
        BOOST_TEST_MESSAGE("In fixture, the seed is " << seed);
        rng.seed(seed);
    }
};

using multitensor::tensor::Tensor;

BOOST_FIXTURE_TEST_SUITE(tests_tensor, fixture)

// Checks the Tensor Class
BOOST_AUTO_TEST_CASE(
    test_tensor,
    *boost::unit_test::tolerance(EPS_PRECISION))
{
    dimension_t nrows = dimension_t(rng() % 10) + 1;
    dimension_t ncols = dimension_t(rng() % 10) + 1;
    dimension_t nlayers = dimension_t(rng() % 10) + 2;

    Tensor<double> T1(nrows, ncols, nlayers);
    BOOST_TEST(T1.size() == nrows * ncols * nlayers);
    BOOST_TEST(T1.order() == 3);

    for (dimension_t i = 0; i < nrows; i++)
    {
        for (dimension_t j = 0; j < ncols; j++)
        {
            for (dimension_t alpha = 0; alpha < nlayers; alpha++)
            {
                BOOST_TEST(T1(i, j, alpha) == 0); // default filled with zero
            }
        }
    }

    dimension_t rand_i = dimension_t(rng() % nrows);
    dimension_t rand_j = dimension_t(rng() % ncols);
    dimension_t rand_alpha = dimension_t(rng() % nlayers);
    double rand_value = double(rng() % 100) - 50;
    T1(rand_i, rand_j, rand_alpha) = rand_value;
    BOOST_TEST(T1(rand_i, rand_j, rand_alpha) == rand_value);

    // Create a tensor with one single layer - basically a matrix
    Tensor<double> T2(nrows, ncols);
    BOOST_TEST(T2.size() == nrows * ncols);
    BOOST_TEST(T2.order() == 2);

    rand_i = dimension_t(rng() % nrows);
    rand_j = dimension_t(rng() % ncols);
    rand_value = double(rng() % 100) - 50;
    T2(rand_i, rand_j) = rand_value;
    BOOST_TEST(T2(rand_i, rand_j) == rand_value);
}

BOOST_AUTO_TEST_SUITE_END()

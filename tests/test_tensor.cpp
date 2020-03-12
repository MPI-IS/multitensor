/*!
 * @file
 *
 * @brief Testing the Tensory class.
 *
 * @author Jean-Claude Passy (jean-claude.passy@tuebingen.mpg.de)
 */

#define BOOST_TEST_MODULE tests_multitensor
#include <random>
#include <ctime>
#include <cstddef>
#include <boost/test/unit_test.hpp>

#include "multitensor/parameters.hpp"
#include "multitensor/tensor.hpp"
#include "fixtures.hpp"

using multitensor::tensor::Tensor;

// Global fixture preparing data
struct fixture_tensor : fixture_rng
{
    dimension_t nrows, ncols, nlayers;
    fixture_tensor()
        : nrows(rng() % 10 + 1),
          ncols(rng() % 10 + 1),
          nlayers(rng() % 10 + 2)
    {
    }
};

BOOST_FIXTURE_TEST_SUITE(tests_tensor, fixture_tensor)

// Checks the Tensor Class
BOOST_AUTO_TEST_CASE(
    test_tensor,
    *boost::unit_test::tolerance(EPS_PRECISION))
{
    // Default tensor with no data
    Tensor<double> T0{};
    BOOST_TEST(T0.size() == 0);

    // Normal tensor
    Tensor<double> T1(nrows, ncols, nlayers);
    BOOST_TEST(T1.size() == nrows * ncols * nlayers);

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

    dimension_t rand_i = rng() % nrows;
    dimension_t rand_j = rng() % ncols;
    dimension_t rand_alpha = rng() % nlayers;
    double rand_value = rng() % 100 - 50;
    T1(rand_i, rand_j, rand_alpha) = rand_value;
    BOOST_TEST(T1(rand_i, rand_j, rand_alpha) == rand_value);

    // Create a tensor with one single layer - basically a matrix
    Tensor<double> T2(nrows, ncols);
    BOOST_TEST(T2.size() == nrows * ncols);

    rand_i = rng() % nrows;
    rand_j = rng() % ncols;
    rand_value = rng() % 100 - 50;
    T2(rand_i, rand_j) = rand_value;
    BOOST_TEST(T2(rand_i, rand_j) == rand_value);
}

// Checks the resizing of a tensor
BOOST_AUTO_TEST_CASE(
    test_resize,
    *boost::unit_test::tolerance(EPS_PRECISION))
{
    // Default tensor with no data
    Tensor<double> T0{};
    BOOST_TEST(T0.size() == 0);

    // Resize
    T0.resize(nrows, ncols, nlayers);
    BOOST_TEST(T0.size() == nrows * ncols * nlayers);

    for (dimension_t i = 0; i < nrows; i++)
    {
        for (dimension_t j = 0; j < ncols; j++)
        {
            for (dimension_t alpha = 0; alpha < nlayers; alpha++)
            {
                BOOST_TEST(T0(i, j, alpha) == 0);
            }
        }
    }

    // Resize again
    dimension_t nrows2 = rng() % 10 + 1;
    dimension_t ncols2 = rng() % 10 + 1;
    T0.resize(nrows2, ncols2);
    BOOST_TEST(T0.size() == nrows2 * ncols2);
}

BOOST_AUTO_TEST_SUITE_END()

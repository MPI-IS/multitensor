// Copyright (c) 2019, Max Planck Society / Software Workshop - Max Planck Institute for Intelligent Systems
// Distributed under the GNU GPL license version 3
// See file LICENSE.md or at https://github.com/MPI-IS/multitensor/LICENSE.md

/*!
 * @file
 *
 * @brief Testing the tensor initializations.
 *
 * @author Jean-Claude Passy (jean-claude.passy@tuebingen.mpg.de)
 */

#include <iostream>
#include <random>
#include <ctime>
#include <vector>
#include <set>
#include <cstdlib>
#include <boost/test/unit_test.hpp>

#include "multitensor/params.hpp"
#include "multitensor/initialization.hpp"
#include "multitensor/tensor.hpp"
#include "fixtures.hpp"

using namespace multitensor::initialization;
using namespace multitensor::tensor;

// Fake generator for random initializations
template <class random_t>
struct fake_generator
{
    const double value;
    fake_generator(random_t &rng)
        : value(double(rng() % 100) - 50)
    {
    }

    double operator()()
    {
        return value;
    }
};

// Global fixture preparing data
struct fixture_init : fixture_rng
{
    const dimension_t nrows, ncols, nlayers;
    fake_generator<std::mt19937_64> fake_rng;
    fixture_init()
        : nrows(rng() % 10 + 1),
          ncols(rng() % 10 + 1),
          nlayers(rng() % 10 + 2),
          fake_rng(rng)
    {
    }
};

BOOST_FIXTURE_TEST_SUITE(tests_init, fixture_init)

// Checks the fully random initialization of a symmetric tensor
BOOST_AUTO_TEST_CASE(test_init_symmetric_tensor_random,
                     *boost::unit_test::tolerance(EPS_PRECISION))
{
    init_symmetric_tensor_random init_obj{};

    // Symmetric tensor
    SymmetricTensor<double> T(nrows, nlayers), Tinit;
    init_obj(Tinit, T, fake_rng);
    BOOST_TEST(T.size() == nrows * nrows * nlayers);
    for (dimension_t alpha = 0; alpha < nlayers; alpha++)
    {
        for (dimension_t j = 0; j < nrows; j++)
        {
            for (dimension_t i = 0; i < nrows; i++)
            {
                BOOST_TEST(T(i, j, alpha) == fake_rng.value);
            }
        }
    }

    // Diagonal tensor
    DiagonalTensor<double> D(nrows, nlayers), Dinit;
    init_obj(Dinit, D, fake_rng);
    BOOST_TEST(D.size() == nrows * nlayers);
    for (dimension_t alpha = 0; alpha < nlayers; alpha++)
    {
        for (dimension_t j = 0; j < nrows; j++)
        {
            BOOST_TEST(D(j, alpha) == fake_rng.value);
        }
    }
}

// Checks the random initialization of symmetric tensor from initial values
BOOST_AUTO_TEST_CASE(test_init_symmetric_tensor_from_initial,
                     *boost::unit_test::tolerance(EPS_PRECISION))
{
    // Symmetric tensor
    SymmetricTensor<double> T(nrows, nlayers), Tinit(nrows, nlayers);
    init_symmetric_tensor_from_initial<SymmetricTensor<double>> init_obj_sym;
    // Initialize values in Tinit
    for (dimension_t alpha = 0; alpha < nlayers; alpha++)
    {
        for (dimension_t j = 0; j < nrows; j++)
        {
            for (dimension_t i = 0; i < nrows; i++)
            {
                Tinit(i, j, alpha) = double(rng() % 100) - 50;
            }
        }
    }
    // Call init
    init_obj_sym(Tinit, T, fake_rng);
    // Checks
    BOOST_TEST(T.size() == nrows * nrows * nlayers);
    for (dimension_t alpha = 0; alpha < nlayers; alpha++)
    {
        for (dimension_t j = 0; j < nrows; j++)
        {
            for (dimension_t i = 0; i < nrows; i++)
            {
                BOOST_TEST(T(i, j, alpha) == (Tinit(i, j, alpha) + EPS_NOISE * fake_rng.value));
            }
        }
    }

    // Diagonal tensor
    DiagonalTensor<double> D(nrows, nlayers), Dinit(nrows, nlayers);
    init_symmetric_tensor_from_initial<DiagonalTensor<double>> init_obj_diag;
    // Initialize values in Tinit
    for (dimension_t alpha = 0; alpha < nlayers; alpha++)
    {
        for (dimension_t j = 0; j < nrows; j++)
        {
            Dinit(j, alpha) = double(rng() % 100) - 50;
        }
    }
    // Call init
    init_obj_diag(Dinit, D, fake_rng);
    // Checks
    BOOST_TEST(D.size() == nrows * nlayers);
    for (dimension_t alpha = 0; alpha < nlayers; alpha++)
    {
        for (dimension_t j = 0; j < nrows; j++)
        {
            BOOST_TEST(D(j, alpha) == (Dinit(j, alpha) + EPS_NOISE * fake_rng.value));
        }
    }
}

// Checks the random initialization of the rows of a matric
BOOST_AUTO_TEST_CASE(test_init_tensor_rows_random,
                     *boost::unit_test::tolerance(EPS_PRECISION))
{
    Matrix<double> M(nrows, ncols);
    const dimension_t nof_elements = rng() % (nrows + 1);
    std::vector<dimension_t> elements(nof_elements);
    for (size_t i = 0; i < nof_elements; i++)
    {
        elements[i] = rng() % nrows;
    }

    // Call init
    init_tensor_rows_random(elements, M, fake_rng);

    // Checks
    std::vector<dimension_t>::iterator it;

    for (dimension_t j = 0; j < ncols; j++)
    {
        for (dimension_t i = 0; i < nrows; i++)
        {
            it = std::find(elements.begin(), elements.end(), i);
            if (it != elements.end())
            {
                BOOST_TEST(M(i, j) == fake_rng.value);
            }
            else
            {
                BOOST_TEST(M(i, j) == 0);
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

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

using multitensor::tensor::DiagonalTensor;
using multitensor::tensor::Matrix;
using multitensor::tensor::SymmetricTensor;
using multitensor::tensor::Tensor;
using multitensor::tensor::Transpose;

// Global fixture preparing data
struct fixture_tensor : fixture_rng
{
    const dimension_t nrows, ncols, nlayers;
    fixture_tensor()
        : nrows(rng() % 10 + 1),
          ncols(rng() % 10 + 1),
          nlayers(rng() % 10 + 2)
    {
    }
};

BOOST_FIXTURE_TEST_SUITE(tests_tensor, fixture_tensor)

// Checks the Tensor class
BOOST_AUTO_TEST_CASE(
    test_tensor,
    *boost::unit_test::tolerance(EPS_PRECISION))
{
    // Default tensor with no data
    Tensor<double> T0{};
    BOOST_TEST(T0.size() == 0);

    // Resize
    T0.resize(nrows, ncols, nlayers);
    BOOST_TEST(T0.size() == nrows * ncols * nlayers);

    // Tensor with zeros
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

    // Tensor from vector
    std::vector<double> vec;
    for (size_t i = 0; i < nrows * ncols * nlayers; i++)
    {
        vec.emplace_back(double(rng() % 100) - 50);
    }
    // Bad dimensions
    BOOST_CHECK_THROW(Tensor<double> T3(nrows + 1, ncols, nlayers, vec),
                      std::runtime_error);
    // OK
    Tensor<double> T3(nrows, ncols, nlayers, vec);
    BOOST_TEST(T3.size() == nrows * ncols * nlayers);
    size_t counter(0);
    for (dimension_t alpha = 0; alpha < nlayers; alpha++)
    {
        for (dimension_t j = 0; j < ncols; j++)
        {
            for (dimension_t i = 0; i < nrows; i++, counter++)
            {
                BOOST_TEST(T3(i, j, alpha) == vec[counter]);
            }
        }
    }

    // Change values
    dimension_t rand_i = rng() % nrows;
    dimension_t rand_j = rng() % ncols;
    dimension_t rand_alpha = rng() % nlayers;
    double rand_value = double(rng() % 100) - 50;
    T1(rand_i, rand_j, rand_alpha) = rand_value;
    BOOST_TEST(T1(rand_i, rand_j, rand_alpha) == rand_value);

    // Get data
    for (size_t i = 0; i < vec.size(); i++)
    {
        BOOST_TEST(vec[i] == T3.get_data()[i]);
    }
}

// Checks the SymmetricTensor class
BOOST_AUTO_TEST_CASE(
    test_symmetric_tensor,
    *boost::unit_test::tolerance(EPS_PRECISION))
{
    // Default symmetric tensor with no data
    SymmetricTensor<double> T0{};
    BOOST_TEST(T0.size() == 0);

    // Resize
    T0.resize(nrows, nlayers);
    BOOST_TEST(T0.size() == nrows * nrows * nlayers);

    // Symmetric tensor with zeros
    SymmetricTensor<double> T1(nrows, nlayers);
    BOOST_TEST(T1.size() == nrows * nrows * nlayers);

    for (dimension_t i = 0; i < nrows; i++)
    {
        for (dimension_t j = 0; j < nrows; j++)
        {
            for (dimension_t alpha = 0; alpha < nlayers; alpha++)
            {
                BOOST_TEST(T1(i, j, alpha) == 0); // default filled with zero
            }
        }
    }

    // Symmetric tensor from vector
    std::vector<double> vec;
    for (size_t i = 0; i < nrows * nrows * nlayers; i++)
    {
        vec.emplace_back(double(rng() % 100) - 50);
    }
    // Bad dimensions
    BOOST_CHECK_THROW(SymmetricTensor<double> T3(nrows + 1, nlayers, vec),
                      std::runtime_error);
    // OK
    SymmetricTensor<double> T3(nrows, nlayers, vec);
    BOOST_TEST(T3.size() == nrows * nrows * nlayers);
    size_t counter(0);
    for (dimension_t alpha = 0; alpha < nlayers; alpha++)
    {
        for (dimension_t j = 0; j < nrows; j++)
        {
            for (dimension_t i = 0; i < nrows; i++, counter++)
            {
                BOOST_TEST(T3(i, j, alpha) == vec[counter]);
            }
        }
    }

    // Change values
    dimension_t rand_i = rng() % nrows;
    dimension_t rand_j = rng() % nrows;
    dimension_t rand_alpha = rng() % nlayers;
    double rand_value = double(rng() % 100) - 50;
    T1(rand_i, rand_j, rand_alpha) = rand_value;
    BOOST_TEST(T1(rand_i, rand_j, rand_alpha) == rand_value);

    // Get data
    for (size_t i = 0; i < vec.size(); i++)
    {
        BOOST_TEST(vec[i] == T3.get_data()[i]);
    }
}

// Checks the Matrix class
BOOST_AUTO_TEST_CASE(
    test_matrix,
    *boost::unit_test::tolerance(EPS_PRECISION))
{
    // Default matrix with no data
    Matrix<double> T0{};
    BOOST_TEST(T0.size() == 0);

    // Resize
    T0.resize(nrows, ncols);
    BOOST_TEST(T0.size() == nrows * ncols);

    // Matrix with zeros
    Matrix<double> T1(nrows, ncols);
    BOOST_TEST(T1.size() == nrows * ncols);

    for (dimension_t i = 0; i < nrows; i++)
    {
        for (dimension_t j = 0; j < ncols; j++)
        {
            BOOST_TEST(T1(i, j) == 0); // default filled with zero
        }
    }

    // Matrix from vector
    std::vector<double> vec;
    for (size_t i = 0; i < nrows * ncols; i++)
    {
        vec.emplace_back(double(rng() % 100) - 50);
    }
    // Bad dimensions
    BOOST_CHECK_THROW(Matrix<double> T3(nrows + 1, ncols, vec),
                      std::runtime_error);
    // OK
    Matrix<double> T3(nrows, ncols, vec);
    BOOST_TEST(T3.size() == nrows * ncols);
    size_t counter(0);
    for (dimension_t j = 0; j < ncols; j++)
    {
        for (dimension_t i = 0; i < nrows; i++, counter++)
        {
            BOOST_TEST(T3(i, j) == vec[counter]);
        }
    }

    // Change values
    dimension_t rand_i = rng() % nrows;
    dimension_t rand_j = rng() % ncols;
    double rand_value = double(rng() % 100) - 50;
    T1(rand_i, rand_j) = rand_value;
    BOOST_TEST(T1(rand_i, rand_j) == rand_value);

    // Get data
    for (size_t i = 0; i < vec.size(); i++)
    {
        BOOST_TEST(vec[i] == T3.get_data()[i]);
    }
}

// Checks the DiagonalTensor class
BOOST_AUTO_TEST_CASE(
    test_diagonal_tensor,
    *boost::unit_test::tolerance(EPS_PRECISION))
{
    // Default diagonal tensor with no data
    DiagonalTensor<double> T0{};
    BOOST_TEST(T0.size() == 0);

    // Resize
    T0.resize(nrows, nlayers);
    BOOST_TEST(T0.size() == nrows * nlayers);

    // Diagonal tensor with zeros
    DiagonalTensor<double> T1(nrows, nlayers);
    BOOST_TEST(T1.size() == nrows * nlayers);

    for (dimension_t i = 0; i < nrows; i++)
    {
        for (dimension_t alpha = 0; alpha < nlayers; alpha++)
        {
            BOOST_TEST(T1(i, alpha) == 0); // default filled with zero
        }
    }

    // Diagonal tensor from vector
    std::vector<double> vec;
    for (size_t i = 0; i < nrows * nlayers; i++)
    {
        vec.emplace_back(double(rng() % 100) - 50);
    }
    // Bad dimensions
    BOOST_CHECK_THROW(DiagonalTensor<double> T3(nrows + 1, nlayers, vec),
                      std::runtime_error);
    // OK
    DiagonalTensor<double> T3(nrows, nlayers, vec);
    BOOST_TEST(T3.size() == nrows * nlayers);
    size_t counter(0);
    for (dimension_t alpha = 0; alpha < nlayers; alpha++)
    {
        for (dimension_t i = 0; i < nrows; i++, counter++)
        {
            BOOST_TEST(T3(i, alpha) == vec[counter]);
        }
    }

    // Change values
    dimension_t rand_i = rng() % nrows;
    dimension_t rand_alpha = rng() % nlayers;
    double rand_value = double(rng() % 100) - 50;
    T1(rand_i, rand_alpha) = rand_value;
    BOOST_TEST(T1(rand_i, rand_alpha) == rand_value);

    // Get data
    for (size_t i = 0; i < vec.size(); i++)
    {
        BOOST_TEST(vec[i] == T3.get_data()[i]);
    }
}

// Checks the Transpose container
BOOST_AUTO_TEST_CASE(
    test_transpose,
    *boost::unit_test::tolerance(EPS_PRECISION))
{
    // Tensor
    Tensor<double> T1(nrows, ncols, nlayers);
    for (dimension_t i = 0; i < nrows; i++)
    {
        for (dimension_t j = 0; j < ncols; j++)
        {
            for (dimension_t alpha = 0; alpha < nlayers; alpha++)
            {
                T1(i, j, alpha) = double(rng() % 100) - 50;
            }
        }
    }
    Transpose T1T(T1);
    BOOST_TEST(std::get<0>(T1.dims()) == std::get<1>(T1T.dims()));
    BOOST_TEST(std::get<1>(T1.dims()) == std::get<0>(T1T.dims()));
    BOOST_TEST(std::get<2>(T1.dims()) == std::get<2>(T1T.dims()));
    for (dimension_t i = 0; i < ncols; i++)
    {
        for (dimension_t j = 0; j < nrows; j++)
        {
            for (dimension_t alpha = 0; alpha < nlayers; alpha++)
            {
                BOOST_TEST(T1T(i, j, alpha) == T1(j, i, alpha));
            }
        }
    }

    // SymmetricTensor
    SymmetricTensor<double> S1(nrows, nlayers);
    for (dimension_t i = 0; i < nrows; i++)
    {
        for (dimension_t j = 0; j < nrows; j++)
        {
            for (dimension_t alpha = 0; alpha < nlayers; alpha++)
            {
                S1(i, j, alpha) = double(rng() % 100) - 50;
            }
        }
    }
    Transpose S1T(S1);
    BOOST_TEST(std::get<0>(S1.dims()) == std::get<1>(S1T.dims()));
    BOOST_TEST(std::get<2>(S1.dims()) == std::get<2>(S1T.dims()));
    for (dimension_t i = 0; i < nrows; i++)
    {
        for (dimension_t j = 0; j < nrows; j++)
        {
            for (dimension_t alpha = 0; alpha < nlayers; alpha++)
            {
                BOOST_TEST(S1T(i, j, alpha) == S1(j, i, alpha));
            }
        }
    }

    // Matrix
    Matrix<double> M1(nrows, ncols);
    for (dimension_t i = 0; i < nrows; i++)
    {
        for (dimension_t j = 0; j < ncols; j++)
        {
            M1(i, j) = double(rng() % 100) - 50;
        }
    }
    Transpose M1T(M1);
    BOOST_TEST(std::get<0>(S1.dims()) == std::get<1>(S1T.dims()));
    BOOST_TEST(std::get<1>(S1.dims()) == std::get<0>(S1T.dims()));
    for (dimension_t i = 0; i < ncols; i++)
    {
        for (dimension_t j = 0; j < nrows; j++)
        {
            BOOST_TEST(M1T(i, j) == M1(j, i));
        }
    }

    // DiagonalTensor
    DiagonalTensor<double> D1(nrows, nlayers);
    for (dimension_t i = 0; i < nrows; i++)
    {
        for (dimension_t alpha = 0; alpha < nlayers; alpha++)
        {
            D1(i, alpha) = double(rng() % 100) - 50;
        }
    }
    Transpose D1T(D1);
    BOOST_TEST(std::get<0>(D1.dims()) == std::get<0>(D1T.dims()));
    BOOST_TEST(std::get<1>(D1.dims()) == std::get<1>(D1T.dims()));
    BOOST_TEST(std::get<2>(D1.dims()) == std::get<2>(D1T.dims()));
    for (dimension_t i = 0; i < nrows; i++)
    {
        for (dimension_t alpha = 0; alpha < nlayers; alpha++)
        {
            BOOST_TEST(D1T(i, alpha) == D1(i, alpha));
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

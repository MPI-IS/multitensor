/*!
 * @file
 *
 * @brief Implementation of initializations functions
 *
 * @author Jean-Claude Passy (jean-claude.passy@tuebingen.mpg.de)
 */

#pragma once

#include "multitensor/parameters.hpp"
#include "multitensor/tensor.hpp"

namespace multitensor
{

//! Initialization functions
namespace initialization
{

/*!
 * @brief Class for the initialization of a symmetric tensor randomly
 */
struct init_tensor_symmetric_random
{
    /*
     * @brief Initializes the tensor
     *
     * @tparam scalar_t Type of the tensor values
     * @tparam random_t Random generator type
     *
     * @param[in,out] T Tensor
     * @param[in,out] random_generator Random generator used for the initialization
     *
     * This function intializes randomly a tensor using the provided random number generator.
     * To facilitate convergence, it is initialized in a symmetric fashion.
     */
    template <class scalar_t,
              class random_t>
    void operator()(tensor::Tensor<scalar_t> &T, random_t &random_generator)
    {
        auto dims = T.dims();
        const dimension_t nrows{std::get<0>(dims)};
        const dimension_t ncols{std::get<1>(dims)};
        const dimension_t ntubes{std::get<2>(dims)};
        assert(nrows == ncols);
        for (size_t alpha = 0; alpha < ntubes; alpha++)
        {
            for (size_t i = 0; i < nrows; i++)
            {
                for (size_t j = i; j < ncols; j++)
                {
                    if (i == j)
                    {
                        T(i, j, alpha) = random_generator();
                    }
                    else
                    {
                        T(i, j, alpha) = T(j, i, alpha) = random_generator();
                    }
                }
            }
        }
    }
};

/*!
 * @brief Class for the partial random initialization of a tensor
 */
struct init_tensor_partial_random
{
    /*
     * @brief Initializes the tensor
     *
     * @tparam scalar_t Type of the tensor values
     * @tparam random_t Random generator type
     *
     * @param[in] Elements Indices of the rows to intialize randomly
     * @param[in,out] T Tensor
     * @param[in,out] random_generator Random generator used for the initialization
     *
     * This function intializes the given rows randomly.
     */
    template <class scalar_t,
              class random_t>
    void operator()(const std::vector<dimension_t> &elements,
                    tensor::Tensor<scalar_t> &T,
                    random_t &random_generator)
    {
        auto dims = T.dims();
        const dimension_t ncols{std::get<1>(dims)};
        const dimension_t ntubes{std::get<2>(dims)};
        for (dimension_t alpha = 0; alpha < ntubes; alpha++)
        {
            for (dimension_t k = 0; k < ncols; k++)
            {
                for (auto j : elements)
                {
                    assert(j < std::get<0>(dims));
                    T(j, k, alpha) = random_generator();
                }
            }
        }
    }
};
} // namespace initialization
} // namespace multitensor
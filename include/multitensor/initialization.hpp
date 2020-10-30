/*!
 * @file
 *
 * @brief Implementation of initializations functions
 *
 * @author Jean-Claude Passy (jean-claude.passy@tuebingen.mpg.de)
 */

#pragma once

#include "multitensor/params.hpp"
#include "multitensor/tensor.hpp"

namespace multitensor
{
//! Initialization functions
namespace initialization
{

//! @brief Class for the fully random initialization of a symmetric tensor
struct init_symmetric_tensor_random
{
    /*
     * @brief Initializes a symmetric tensor randomly
     *
     * @tparam tensor_t Tensor type
     * @tparam random_t Random generator type
     *
     * @param[in] Tinit Tensor used to store initial values (not used here)
     * @param[in,out] T Tensor
     * @param[in,out] random_generator Random generator used for the initialization
     *
     * This function intializes randomly a tensor using the provided random number generator.
     * Initialize both
     *     - a @ref Tensor (initialized in a symmetric fashion to facilitate convergencea)
     *     - a @ref DiagonalTensor.
     *
     * @pre Tensor must have been resized appropriately.
     */
    template <class tensor_t,
              class random_t>
    void operator()(const tensor_t & /*Tinit*/, tensor_t &T, random_t &random_generator)
    {
        using scalar_t = std::decay_t<decltype(T(0, 0))>;

        auto dims = T.dims();
        const dimension_t nrows{std::get<0>(dims)};
        const dimension_t ncols{std::get<1>(dims)};
        const dimension_t ntubes{std::get<2>(dims)};
        for (size_t alpha = 0; alpha < ntubes; alpha++)
        {
            for (size_t i = 0; i < nrows; i++)
            {
                // Compile-time if statement
                // DiagonalTensor
                if constexpr (std::is_same_v<tensor_t, tensor::DiagonalTensor<double>>)
                {
                    T(i, alpha) = static_cast<scalar_t>(random_generator());
                }
                else // SymmetricTensor
                {
                    for (size_t j = i; j < ncols; j++)
                    {
                        if (i == j)
                        {
                            T(i, j, alpha) = static_cast<scalar_t>(random_generator());
                        }
                        else
                        {
                            T(i, j, alpha) = T(j, i, alpha) = static_cast<scalar_t>(random_generator());
                        }
                    }
                }
            }
        }
    }
};

/*
 * @brief Class for initializing of a symmetric tensor from initial values
 *
 * @tparam tensor_t Tensor type
 */
template <class tensor_t>
struct init_symmetric_tensor_from_initial
{
    //! Tensor storing the initial values that will be used for each realizations
    tensor_t tensor_init{};

    /*
     * @brief Initializes a symmetric tensor from initial values
     *
     * @tparam random_t Random generator type
     *
     * @param[in,out] Tinit Tensor used to store initial values
     * @param[in,out] T Tensor
     * @param[in,out] random_generator Random generator used for the initialization
     *
     * This function intializes the tensor from intial values and add some noise afterwards.
     * Initialize both
     *     - a @ref Tensor (initialized in a symmetric fashion to facilitate convergencea)
     *     - a @ref DiagonalTensor.
     *
     * @pre Tensor must have been resized appropriately.
     */
    template <class random_t>
    void operator()(const tensor_t &Tinit, tensor_t &T, random_t &random_generator)
    {
        const dimension_t nof_rows(std::get<0>(Tinit.dims()));
        const dimension_t nof_tubes(std::get<2>(Tinit.dims()));

        // If it has not been done, save initial values
        if (tensor_init.size() == 0)
        {
            tensor_init = Tinit;
        }

        // Copy initial values
        T = tensor_init;

        // Add noise
        for (size_t alpha = 0; alpha < nof_tubes; alpha++)
        {
            for (size_t k = 0; k < nof_rows; k++)
            {
                // DiagonalTensor
                if constexpr (std::is_same_v<tensor_t, tensor::DiagonalTensor<double>>)
                {
                    T(k, alpha) += EPS_NOISE * random_generator();
                }
                else // SymmetricTensor
                {
                    for (size_t q = 0; q < nof_rows; q++)
                    {
                        T(k, q, alpha) += EPS_NOISE * random_generator();
                    }
                }
            }
        }
    }
};

/*
 * @brief Initializes the rows of a matrix randomly
 *
 * @tparam scalar_t Type of the tensor values
 * @tparam random_t Random generator type
 *
 * @param[in] Elements Indices of the rows to intialize randomly
 * @param[in,out] mat Matrix
 * @param[in,out] random_generator Random generator used for the initialization
 *
 * This function intializes the given rows randomly.
 *
 * @pre Matrix must have been resized appropriately.
 */
template <class scalar_t,
          class random_t>
void init_tensor_rows_random(const std::vector<dimension_t> &elements,
                             tensor::Matrix<scalar_t> &mat,
                             random_t &random_generator)
{
    auto dims = mat.dims();
    const dimension_t ncols{std::get<1>(dims)};
    for (dimension_t k = 0; k < ncols; k++)
    {
        for (auto j : elements)
        {
            assert(j < std::get<0>(dims));
            mat(j, k) = static_cast<scalar_t>(random_generator());
        }
    }
}

} // namespace initialization
} // namespace multitensor

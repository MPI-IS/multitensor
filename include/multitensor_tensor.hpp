/*!
 * @file
 *
 * @brief Classes and functions used for algebric manipulations.
 *
 * @author Jean-Claude Passy (jean-claude.passy@tuebingen.mpg.de)
 */

#pragma once

#include <vector>

#include "multitensor_parameters.hpp"

namespace multitensor
{
namespace tensor
{

/*!
 * @brief Third-order tensor
 *
 * Class for a tensor of order 3.
 */
template <class scalar_t>
class Tensor
{
private:
    dimension_t nrows;
    dimension_t ncols;
    dimension_t ntubes; // 3rd component
    std::vector<scalar_t> data;

public:
    /*!
     * @brief Constructor of a tensor
     *
     * @param[in] nrows_ number of rows
     * @param[in] ncols_ number of columns
     * @param[in] ntubes_ number of tubes/layers
     */
    Tensor(dimension_t nrows_, dimension_t ncols_, dimension_t ntubes_ = 1)
        : nrows(nrows_),
          ncols(ncols_),
          ntubes(ntubes_),
          data(nrows * ncols * ntubes, scalar_t(0))
    {
    }

    /*!
     * @brief Initialize a symmetric tensor randomly
     *
     *
     */
    template <class random_t>
    void randomize_symmetric(random_t &random_generator)
    {
        assert(nrows == ncols);
        for (size_t alpha = 0; alpha < ntubes; alpha++)
        {
            for (size_t i = 0; i < nrows; i++)
            {
                for (size_t j = i; j < ncols; j++)
                {
                    if (i == j)
                    {
                        operator()(i, j, alpha) = random_generator();
                    }
                    else
                    {
                        operator()(i, j, alpha) = operator()(j, i, alpha) = random_generator();
                    }
                }
            }
        }
    }

    /*!
     * @brief Initialize specific values of a tensor randomly
     *
     *
     */
    template <class random_t>
    void randomize_partial(random_t &random_generator, const std::vector<size_t> &elements)
    {
        for (size_t alpha = 0; alpha < ntubes; alpha++)
        {
            for (size_t k = 0; k < ncols; k++)
            {
                /*
                for (auto j : elements)
                {
                    operator()(j, k, alpha) = random_generator();
                }
                */
                for (size_t i = 0; i < elements.size(); i++)
                {
                    operator()(elements[i], k, alpha) = random_generator();
                }
            }
        }
    }

    /*!
     * @brief Returns a reference to an element of the tensor
     *
     * @param[in] i row index
     * @param[in] j column index
     * @param[in] alpha tube/layer index
     */
    scalar_t &operator()(const size_t i, const size_t j, const size_t alpha = 0)
    {
        return data[get_index(i, j, alpha)];
    }

    /*!
     * @brief Returns a const reference to an element in the tensor
     *
     * @param[in] i row index
     * @param[in] j column index
     * @param[in] alpha tube/layer index
     */
    const scalar_t &operator()(const size_t i, const size_t j, const size_t alpha = 0) const
    {
        return data[get_index(i, j, alpha)];
    }

    //! @brief Returns the size of the tensor.
    size_t size() const noexcept
    {
        return data.size();
    }

    //! @brief Returns the dimensions of a tensor
    auto dims() const noexcept
    {
        return std::make_tuple(nrows, ncols, ntubes);
    }

    //! @brief Returns the order of the tensor.
    size_t order() const noexcept
    {
        return (ntubes > 1) ? 3
                            : 2;
    }

private:
    //! @brief Helper returning the index value corresponding to a tensor element
    size_t get_index(const size_t i, const size_t j, const size_t alpha) const
    {
        assert(i < nrows);
        assert(j < ncols);
        assert(alpha < ntubes);

        size_t index = alpha * ncols * nrows + j * nrows + i;
        assert(index < size());
        return index;
    }
};

} // namespace tensor
} // namespace multitensor

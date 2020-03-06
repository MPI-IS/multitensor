/*!
 * @file
 *
 * @brief Classes and functions used for algebric manipulations.
 *
 * @author Jean-Claude Passy (jean-claude.passy@tuebingen.mpg.de)
 */

#pragma once

#include <vector>
#include <fstream>
#include <iomanip>
#include <cstddef>
#include <boost/filesystem.hpp>

#include "multitensor/parameters.hpp"

namespace multitensor
{

//! Tensor and algebric manipulations
namespace tensor
{

/*!
 * @brief Class for a tensor of order 3.
 *
 * @tparam scalar_t Type of the tensor values
 */
template <class scalar_t>
class Tensor
{
private:
    dimension_t nrows;
    dimension_t ncols;
    dimension_t ntubes; // 3rd component
    std::vector<scalar_t> data;

    /*!
     * @brief Index value corresponding to a tensor element
     *
     * @param[in] i Row index
     * @param[in] j Column index
     * @param[in] alpha Tube/layer index
     *
     * @returns Index number
     */
    size_t get_index(const size_t i, const size_t j, const size_t alpha) const
    {
        assert(i < nrows);
        assert(j < ncols);
        assert(alpha < ntubes);

        size_t index = alpha * ncols * nrows + j * nrows + i;
        assert(index < size());
        return index;
    }

public:
    /*!
     * @brief Tensor constructor
     *
     * @param[in] nrows Number of rows
     * @param[in] ncols Number of columns
     * @param[in] ntubes Number of tubes/layers
     */
    Tensor(dimension_t nrows, dimension_t ncols, dimension_t ntubes = 1)
        : nrows(nrows),
          ncols(ncols),
          ntubes(ntubes),
          data(nrows * ncols * ntubes, scalar_t(0))
    {
    }

    /*!
     * @brief Initialize a symmetric tensor randomly
     *
     * @tparam random_t Random generator type
     *
     * @param[in,out] random_generator Random generator used for the initialization
     *
     * This function intializes randomly a tensor using the provided random number generator.
     * To facilitate convergence, it is initialized in a symmetric fashion.
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
     * @tparam random_t Random generator type
     *
     * @param[in] elements Indices of elements to initialize randomly
     * @param[in,out] random_generator Random generator used for the initialization
     */
    template <class random_t>
    void randomize_partial(const std::vector<size_t> &elements, random_t &random_generator)
    {
        for (size_t alpha = 0; alpha < ntubes; alpha++)
        {
            for (size_t k = 0; k < ncols; k++)
            {
                for (auto j : elements)
                {
                    operator()(j, k, alpha) = random_generator();
                }
            }
        }
    }

    /*!
     * @brief Returns a reference to an element of the tensor
     *
     * @param[in] i Row index
     * @param[in] j Column index
     * @param[in] alpha Tube/layer index
     *
     * @returns Tensor element
     */
    scalar_t &operator()(const size_t i, const size_t j, const size_t alpha = 0)
    {
        return data[get_index(i, j, alpha)];
    }

    /*!
     * @brief Returns a const reference to an element in the tensor
     *
     * @param[in] i Row index
     * @param[in] j Column index
     * @param[in] alpha Tube/layer index
     *
     * @returns Tensor element
     */
    const scalar_t &operator()(const size_t i, const size_t j, const size_t alpha = 0) const
    {
        return data[get_index(i, j, alpha)];
    }

    //! @brief Returns the size of the tensor
    size_t size() const noexcept
    {
        return data.size();
    }

    /*!
     * @brief The dimensions of a tensor
     *
     * @returns Tuple with the three dimensions
     */
    auto dims() const noexcept
    {
        return std::make_tuple(nrows, ncols, ntubes);
    }
};

} // namespace tensor
} // namespace multitensor
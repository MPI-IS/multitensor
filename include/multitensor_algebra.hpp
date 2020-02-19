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
namespace algebra
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
    dimension_t nlayers; // 3rd component, also called tube
    std::vector<scalar_t> data;

public:
    /*! Constructor of a tensor
     *
     * @param[in] nrows_ number of rows
     * @param[in] ncols_ number of columns
     * @param[in] nlayers_ number of layers
     */
    Tensor(dimension_t nrows_, dimension_t ncols_, dimension_t nlayers_ = 1)
        : nrows(nrows_),
          ncols(ncols_),
          nlayers(nlayers_),
          data(nrows * ncols * nlayers, scalar_t(0))
    {
    }

    /*! Returns a reference to an element of the tensor
     *
     * @param[in] i row index
     * @param[in] j column index
     * @param[in] alpha layer index
     */
    scalar_t &operator()(const size_t i, const size_t j, const size_t alpha = 0)
    {
        return data[get_index(i, j, alpha)];
    }

    /*! Returns a const reference to an element in the tensor
     *
     * @param[in] i row index
     * @param[in] j column index
     * @param[in] alpha layer index
     */
    const scalar_t &operator()(const size_t i, const size_t j, const size_t alpha = 0) const
    {
        return data[get_index(i, j, alpha)];
    }

    //! @brief The size of the tensor.
    size_t size() const noexcept { return this->data.size(); }

private:
    // Helper returning the index value corresponding to a tensor element
    size_t get_index(const size_t i, const size_t j, const size_t alpha) const
    {
        assert(i < nrows);
        assert(j < ncols);
        assert(alpha < nlayers);

        size_t index = alpha * ncols * nrows + i * ncols + j;
        assert(index < size());
        return index;
    }
};

} // namespace algebra
} // namespace multitensor

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
#include <cassert>

#include "multitensor/parameters.hpp"

namespace multitensor
{

//! Tensor and algebric manipulations
namespace tensor
{

/*!
 * @brief Base class for a tensor of order 3.
 *
 * @tparam scalar_t Type of the tensor values
 */
template <class scalar_t>
class Tensor
{
protected:
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
    //! @brief Default constructor
    Tensor() = default;

    /*!
     * @brief Tensor constructor
     *
     * @param[in] nrows Number of rows
     * @param[in] ncols Number of columns
     * @param[in] ntubes Number of tubes/layers
     */
    Tensor(dimension_t nrows, dimension_t ncols, dimension_t ntubes)
        : Tensor()
    {
        resize(nrows, ncols, ntubes);
    }

    /*!
     * @brief Tensor constructor from a vector
     *
     * @param[in] nrows Number of rows
     * @param[in] ncols Number of columns
     * @param[in] ntubes Number of tubes/layers
     * @param[in] vec Data
     */
    Tensor(dimension_t nrows, dimension_t ncols, dimension_t ntubes, const std::vector<scalar_t> &vec)
        : nrows(nrows),
          ncols(ncols),
          ntubes(ntubes),
          data(vec)
    {

        if (data.size() != nrows * ncols * ntubes)
        {
            throw std::runtime_error(
                "[multitensor::tensor] Dimensions inconsistent with vector size: " +
                std::to_string(nrows) + "x" + std::to_string(ncols) + "x" +
                std::to_string(ntubes) + " != " + std::to_string(data.size()) + "\n");
        }
    }

    /*!
     * @brief Function resizing the tensor
     *
     * @param[in] nrows Number of rows
     * @param[in] ncols Number of columns
     * @param[in] ntubes Number of tubes/layers
     */
    void resize(dimension_t nrows_, dimension_t ncols_, dimension_t ntubes_)
    {
        nrows = nrows_;
        ncols = ncols_;
        ntubes = ntubes_;

        data.clear();
        data.assign(nrows * ncols * ntubes, scalar_t(0));
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

    //! @brief Access the data vector
    const std::vector<scalar_t> &get_data() const noexcept
    {
        return data;
    }
}; // end class tensor

/*!
 * @brief Class for a symmetric tensor.
 *
 * @tparam scalar_t Type of the tensor values
 */
template <class scalar_t>
class SymmetricTensor : public Tensor<scalar_t>
{
public:
    //! @brief Default constructor
    using Tensor<scalar_t>::Tensor;

    /*!
     * @brief SymmetricTensor constructor
     *
     * @param[in] nrows Number of rows/columns
     * @param[in] ntubes Number of tubes/layers
     */
    SymmetricTensor(dimension_t nrows, dimension_t ntubes)
        : Tensor<scalar_t>::Tensor(nrows, nrows, ntubes)
    {
    }

    /*!
     * @brief SymmetricTensor constructor from a vector
     *
     * @param[in] nrows Number of rows/columns
     * @param[in] ntubes Number of tubes/layers
     * @param[in] vec Data
     */
    SymmetricTensor(dimension_t nrows, dimension_t ntubes, const std::vector<scalar_t> &vec)
        : Tensor<scalar_t>::Tensor(nrows, nrows, ntubes, vec)
    {
    }

    /*!
     * @brief Function resizing the tensor
     *
     * @param[in] nrows Number of rows/columns
     * @param[in] ntubes Number of tubes/layers
     */
    void resize(dimension_t nrows, dimension_t ntubes)
    {
        Tensor<scalar_t>::resize(nrows, nrows, ntubes);
    }
};

/*!
 * @brief Class for a matrix
 *
 * @tparam scalar_t Type of the tensor values
 *
 * @note This is a single layer tensor
 */
template <class scalar_t>
class Matrix : public Tensor<scalar_t>
{
public:
    //! @brief Default constructor
    using Tensor<scalar_t>::Tensor;

    /*!
     * @brief Matrix constructor
     *
     * @param[in] nrows Number of rows
     * @param[in] ncols Number of columns
     */
    Matrix(dimension_t nrows, dimension_t ncols)
        : Tensor<scalar_t>::Tensor(nrows, ncols, 1)
    {
    }

    /*!
     * @brief Matrix constructor from a vector
     *
     * @param[in] nrows Number of rows
     * @param[in] ncols Number of columns
     * @param[in] vec Data
     */
    Matrix(dimension_t nrows, dimension_t ncols, const std::vector<scalar_t> &vec)
        : Tensor<scalar_t>::Tensor(nrows, ncols, 1, vec)
    {
    }

    /*!
     * @brief Function resizing the matrix
     *
     * @param[in] nrows Number of rows/columnsrows
     * @param[in] ncols Number of columns
     */
    void resize(dimension_t nrows, dimension_t ncols)
    {
        Tensor<scalar_t>::resize(nrows, ncols, 1);
    }

    /*!
     * @brief Returns a reference to an element of the matrix
     *
     * @param[in] i Row index
     * @param[in] j Column index
     *
     * @returns Tensor element
     */
    scalar_t &operator()(const size_t i, const size_t j)
    {
        return Tensor<scalar_t>::operator()(i, j, 0);
    }

    /*!
     * @brief Returns a const reference to an element in the matrix
     *
     * @param[in] i Row index
     * @param[in] j Column index
     *
     * @returns Tensor element
     */
    const scalar_t &operator()(const size_t i, const size_t j) const
    {
        return Tensor<scalar_t>::operator()(i, j, 0);
    }
};

/*!
 * @brief Class for a diagonal tensor
 *
 * @tparam scalar_t Type of the tensor values
 *
 * @note The values are non-0 only on the diagonals
 */
template <class scalar_t>
class DiagonalTensor : public Tensor<scalar_t>
{
public:
    //! @brief Default constructor
    using Tensor<scalar_t>::Tensor;

    /*!
     * @brief DiagonalTensor constructor
     *
     * @param[in] nrows Number of rows/columns
     * @param[in] ntubes Number of tubes/layers
     */
    DiagonalTensor(dimension_t nrows, dimension_t ntubes)
        : Tensor<scalar_t>::Tensor(nrows, 1, ntubes)
    {
    }

    /*!
     * @brief DiagonalTensor constructor from a vector
     *
     * @param[in] nrows Number of rows/columns
     * @param[in] ntubes Number of tubes/layers
     * @param[in] vec Data
     */
    DiagonalTensor(dimension_t nrows, dimension_t ntubes, const std::vector<scalar_t> &vec)
        : Tensor<scalar_t>::Tensor(nrows, 1, ntubes, vec)
    {
    }

    /*!
     * @brief Function resizing the DiagonalTensor
     *
     * @param[in] nrows Number of rows/columns
     * @param[in] ntubes Number of tubes/layers
     */
    void resize(dimension_t nrows, dimension_t ntubes)
    {
        Tensor<scalar_t>::resize(nrows, 1, ntubes);
    }

    /*!
     * @brief Returns a reference to an element of the DiagonalTensor
     *
     * @param[in] i Row index
     * @param[in] alpha Tube/layer index
     *
     * @returns Tensor element
     */
    scalar_t &operator()(const size_t i, const size_t alpha)
    {
        return Tensor<scalar_t>::operator()(i, 0, alpha);
    }

    /*!
     * @brief Returns a const reference to an element in the DiagonalTensor
     *
     * @param[in] i Row index
     * @param[in] alpha Tube/layer index
     *
     * @returns Tensor element
     */
    const scalar_t &operator()(const size_t i, const size_t alpha) const
    {
        return Tensor<scalar_t>::operator()(i, 0, alpha);
    }
};

/*!
 * @brief Container for transposing a tensor
 *
 * @tparam tensor_t Tensor type
 */
template <class tensor_t>
class Transpose
{
private:
    tensor_t &tensor;

public:
    /*
     * @brief Constructor
     *
     * param[in] tensor Tensor to transpose
     */
    Transpose(tensor_t &tensor)
        : tensor(tensor)
    {
    }

    using scalar_t = std::decay_t<decltype(tensor(0, 0))>;

    /*!
     * @brief Returns a reference to an element of the transposed tensor
     *
     * @param[in] i Row index
     * @param[in] j Column index
     * @param[in] alpha Tube/layer index
     *
     * @returns Tensor element
     */
    scalar_t &operator()(const size_t i, const size_t j, const size_t alpha)
    {
        return tensor(j, i, alpha);
    }
    /*!
     * @brief Returns a const reference to an element in the transposed tensor
     *
     * @param[in] i Row index
     * @param[in] j Column index
     * @param[in] alpha Tube/layer index
     *
     * @returns Tensor element
     */
    const scalar_t &operator()(const size_t i, const size_t j, const size_t alpha) const
    {
        return tensor(j, i, alpha);
    }

    /*!
     * @brief Returns a reference to an element of the transposed tensor
     *
     * @param[in] i Row index
     * @param[in] j Column index
     *
     * @returns Tensor element
     */
    scalar_t &operator()(const size_t i, const size_t j)
    {
        // Tr(DiagTensor) = DiagTensor
        if constexpr (std::is_same_v<tensor_t, DiagonalTensor<double>>)
        {
            return tensor(i, j);
        }
        return tensor(j, i);
    }
    /*!
     * @brief Returns a const reference to an element in the transposed tensor
     *
     * @param[in] i Row index
     * @param[in] j Column index
     *
     * @returns Tensor element
     */
    const scalar_t &operator()(const size_t i, const size_t j) const
    {
        // Tr(DiagTensor) = DiagTensor
        if constexpr (std::is_same_v<tensor_t, DiagonalTensor<double>>)
        {
            return tensor(i, j);
        }
        return tensor(j, i);
    }

    /*!
     * @brief The dimensions of the transpose
     *
     * @returns Tuple with the three dimensions
     */
    auto dims() const noexcept
    {
        // Tr(DiagTensor) = DiagTensor
        if constexpr (std::is_same_v<tensor_t, DiagonalTensor<double>>)
        {
            return tensor.dims();
        }
        return std::make_tuple(std::get<1>(tensor.dims()),
                               std::get<0>(tensor.dims()),
                               std::get<2>(tensor.dims()));
    }
};

} // namespace tensor
} // namespace multitensor
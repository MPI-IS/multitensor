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

#include "multitensor/parameters.hpp"

namespace multitensor
{

//! Tensor and algebric manipulations
namespace tensor
{

/*!
 * @brief Class for a tensor of order 3.
 *
 * @tparam scalar_t the type of the tensor values
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
     * @param[in] nrows_ Number of rows
     * @param[in] ncols_ Number of columns
     * @param[in] ntubes_ Number of tubes/layers
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

    /*!
     * @brief Write tensor into a file
     *
     * @param[in] filename Name of the output file
     */
    void write_affinity_file(const std::string &output_filename,
                             const double &L2,
                             const unsigned int &nof_realizations)
    {
        std::ofstream stream_out(output_filename.c_str());
        if (stream_out.fail())
        {
            throw(std::runtime_error(
                std::string("In write_tensor_file, failed to open ") + output_filename));
        }

        // Likelihood and number of realizations
        stream_out << "# Max likelihood= " << static_cast<int>(L2)
                   << " N_real=" << nof_realizations << std::endl;

        stream_out << std::setprecision(6);
        for (size_t alpha = 0; alpha < ntubes; alpha++)
        {
            stream_out << "a= " << alpha << std::endl;
            for (size_t k = 0; k < nrows; k++)
            {
                for (size_t q = 0; q < ncols; q++)
                {
                    stream_out << operator()(k, q, alpha) << " ";
                }
                stream_out << std::endl;
            }
            stream_out << std::endl;
        }
    }

    /*!
     * @brief Write tensor into a file
     *
     * @param[in] filename Name of the output file
     *
     * @todo: Improved using report
     */
    template <class network_t>
    void write_membership_file(const std::string &output_filename,
                               const network_t &A,
                               const double &L2,
                               const unsigned int &nof_realizations)
    {
        std::ofstream stream_out(output_filename.c_str());
        if (stream_out.fail())
        {
            throw(std::runtime_error(
                std::string("In write_tensor_file, failed to open ") + output_filename));
        }

        // Likelihood and number of realizations
        stream_out << "# Max likelihood= " << static_cast<int>(L2)
                   << " N_real=" << nof_realizations << std::endl;

        stream_out << std::setprecision(6);
        for (size_t alpha = 0; alpha < ntubes; alpha++)
        {
            for (size_t k = 0; k < nrows; k++)
            {
                stream_out << A(0)[k].label << " ";
                for (size_t q = 0; q < ncols; q++)
                {
                    stream_out << operator()(k, q, alpha) << " ";
                }
                stream_out << std::endl;
            }
        }
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

    /*!
     * @brief The order of the tensor
     *
     * @returns 2 if the number of layers == 1, otherwise 3
     */
    size_t order() const noexcept
    {
        return (ntubes > 1) ? 3 : 2;
    }
};

} // namespace tensor
} // namespace multitensor
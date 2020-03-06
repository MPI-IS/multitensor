/*!
 * @file
 *
 * @brief Utility classes and functions.
 *
 * @author Jean-Claude Passy (jean-claude.passy@tuebingen.mpg.de)
 */

#pragma once

#include <set>

#include "multitensor/tensor.hpp"

namespace multitensor
{
//! Utilities
namespace utils
{
//! @brief Class for containing the results of the algorithm
struct Report
{
    //! The number of realizations
    unsigned int nof_realizations;

    //! The likelihood of the best results
    std::vector<double> vec_L2;

    //! The number of iterations for each realization
    std::vector<size_t> vec_iter;

    //! The reason for termination
    std::vector<const char *> vec_term_reason;

    //! Duration (in seconds) of the full run
    double duration;

    Report()
        : nof_realizations(0),
          vec_L2(0),
          vec_iter(0),
          vec_term_reason(0),
          duration(0)
    {
    }

    double max_L2() const noexcept
    {
        if (vec_L2.size() == 0)
        {
            return std::numeric_limits<double>::lowest();
        }
        else
        {
            auto i = std::max_element(vec_L2.begin(), vec_L2.end());
            return vec_L2[std::distance(vec_L2.begin(), i)];
        }
    }
};

/*!
 * @brief Utility function calculating the number of vertices from input data
 *
 * @param[in] edges_start Labels of vertices where an edge starts
 * @param[in] edges_end Labels of vertices where an edge ends
 *
 * @returns Number of vertices
 */
template <class vertex_t>
size_t calculate_num_vertices(const std::vector<vertex_t> &edges_start,
                              const std::vector<vertex_t> &edges_end)
{
    std::set<vertex_t> set_e_in(edges_start.begin(), edges_start.end());
    std::set<vertex_t> set_e_out(edges_end.begin(), edges_end.end());
    std::set<vertex_t> set_e;
    std::merge(set_e_in.begin(), set_e_in.end(),
               set_e_out.begin(), set_e_out.end(),
               std::inserter(set_e, set_e.begin()));
    return set_e.size();
}

/*!
 * @brief Write affinity file
 *
 * @param[in] output_filename Name of the output file
 * @param[in] w Affinity tensor
 * @param[in] results Results of the run
 */
template <class scalar_t>
void write_affinity_file(const boost::filesystem::path &output_filename,
                         const tensor::Tensor<scalar_t> &w,
                         const utils::Report &results)
{
    std::ofstream stream_out(output_filename.string());
    if (stream_out.fail())
    {
        throw(std::runtime_error(
            std::string("In write_affinity_file, failed to open ") + output_filename.string()));
    }

    // Likelihood and number of realizations
    stream_out << "# Max likelihood= " << static_cast<int>(results.max_L2())
               << " N_real=" << results.nof_realizations << std::endl;

    stream_out << std::setprecision(6);
    auto dims = w.dims();
    const dimension_t nrows{std::get<0>(dims)};
    const dimension_t ncols{std::get<1>(dims)};
    const dimension_t ntubes{std::get<2>(dims)};
    for (dimension_t alpha = 0; alpha < ntubes; alpha++)
    {
        stream_out << "a= " << alpha << std::endl;
        for (dimension_t k = 0; k < nrows; k++)
        {
            for (dimension_t q = 0; q < ncols; q++)
            {
                stream_out << w(k, q, alpha) << " ";
            }
            stream_out << std::endl;
        }
        stream_out << std::endl;
    }
}

/*!
 * @brief Write tensor into a file
 *
 * @param[in] output_filename Name of the output file
 * @param[in] A Network
 * @param[in] t Tensor (u or v)
 * @param[in] results Results of the run
 */
template <class scalar_t,
          class network_t>
void write_membership_file(const boost::filesystem::path &output_filename,
                           const network_t &A,
                           const tensor::Tensor<scalar_t> &t,
                           const utils::Report &results)
{
    std::ofstream stream_out(output_filename.string());
    if (stream_out.fail())
    {
        throw(std::runtime_error(
            std::string("In write_membership_file, failed to open ") + output_filename.string()));
    }

    // Likelihood and number of realizations
    stream_out << "# Max likelihood= " << static_cast<int>(results.max_L2())
               << " N_real=" << results.nof_realizations << std::endl;

    stream_out << std::setprecision(6);
    auto dims = t.dims();
    const dimension_t nrows{std::get<0>(dims)};
    const dimension_t ncols{std::get<1>(dims)};
    const dimension_t ntubes{std::get<2>(dims)};
    for (size_t alpha = 0; alpha < ntubes; alpha++)
    {
        for (size_t k = 0; k < nrows; k++)
        {
            stream_out << A(0)[k].label << " ";
            for (size_t q = 0; q < ncols; q++)
            {
                stream_out << t(k, q, alpha) << " ";
            }
            stream_out << std::endl;
        }
    }
}

/*!
 * @brief Write info file
 *
 * @param[in] output_filename Name of the output file
 * @param[in] results Results of the run
 */
inline void write_info_file(const boost::filesystem::path &output_filename,
                            const Report &results)
{
    std::ofstream stream_out(output_filename.string());
    if (stream_out.fail())
    {
        throw(std::runtime_error(
            std::string("In write_info_file, failed to open ") + output_filename.string()));
    }

    stream_out << "# Number of realization = " << results.nof_realizations << std::endl;
    stream_out << "# Maximum Likelihood = " << results.max_L2() << std::endl;
    stream_out << "# Duration (s) = " << results.duration << std::endl;
    stream_out << "# real  num_iters  term_reason  L2" << std::endl;
    for (size_t i = 0; i < results.vec_iter.size(); i++)
    {
        stream_out << "  " << i << "  " << results.vec_iter[i] << "  " << results.vec_term_reason[i]
                   << "  " << results.vec_L2[i] << std::endl;
    }
}
} // namespace utils
} // namespace multitensor

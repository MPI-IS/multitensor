/*!
 * @file
 *
 * @brief Application utilities.
 *
 * @author Jean-Claude Passy (jean-claude.passy@tuebingen.mpg.de)
 */

#pragma once

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cassert>
#include <algorithm>
#include <cstdlib>
#include <functional>
#include <numeric>
#include <string>
#include <vector>
#include <set>
#include <boost/filesystem.hpp>

#include "multitensor/parameters.hpp"
#include "multitensor/tensor.hpp"
#include "multitensor/utils.hpp"

using namespace multitensor;

//! @brief Interpret a command option given by the user
char *
get_cmd_option(char **begin, char **end, const std::string &option);

//! @brief Check if a command option has been specified by the user
bool cmd_option_exists(char **begin, char **end, const std::string &option);

/*!
 * @brief Read adjacency matrix data and creates the necessary vectors
 *
 * @param[in] filename Name of the file containing the data
 * @param[in,out] edges_start Labels of vertices where an edge starts
 * @param[in,out] edges_end Labels of vertices where an edge ends
 * @param[in,out] edges_weight Edge weights
 */
template <class scalar_t>
void read_adjacency_data(const std::string &filename,
                         std::vector<size_t> &edges_start,
                         std::vector<size_t> &edges_end,
                         std::vector<scalar_t> &edges_weight)
{
    edges_start.clear();
    edges_end.clear();
    edges_weight.clear();
    std::ifstream in(filename.c_str());
    if (in.fail())
    {
        throw std::runtime_error(
            std::string("In read_adjacency_data, failed to open ") + filename);
    }

    std::cout << "Reading adjacency file " << filename << std::endl;

    std::string line;
    while (!in.eof())
    {
        std::getline(in, line);
        if (line.size() == 0)
            continue; // skip over empty lines

        // Remove trailing whitespaces
        line.erase(line.find_last_not_of(" ") + 1);

        std::vector<scalar_t> current_weights;
        std::istringstream is(line);

        // First character should be an E for edge
        std::string tok;
        is >> tok;
        assert(tok == "E");

        // Second two characters are the edges ids
        // TO DO: might use tuple later on
        size_t current_edge_in, current_edge_out;
        is >> current_edge_in >> current_edge_out;

        // Read the rest of the data
        scalar_t value;
        while (is >> value)
        {

            current_weights.push_back(value);
        }

        // Insert data at the end of the vectors
        edges_start.push_back(current_edge_in);
        edges_end.push_back(current_edge_out);
        edges_weight.insert(
            std::end(edges_weight), std::begin(current_weights), std::end(current_weights));
    }

    // Close file
    in.close();
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
 * @brief Write info file
 *
 * @param[in] output_filename Name of the output file
 * @param[in] results Results of the run
 */
void write_info_file(const boost::filesystem::path &output_filename,
                     const utils::Report &results);

/*!
 * @brief Write tensor into a file
 *
 * @param[in] output_filename Name of the output file
 * @param[in] labels Vertices labels
 * @param[in] t Tensor (u or v)
 * @param[in] results Results of the run
 */
template <class scalar_t>
void write_membership_file(const boost::filesystem::path &output_filename,
                           const std::vector<size_t> &labels,
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
            stream_out << labels[k] << " ";
            for (size_t q = 0; q < ncols; q++)
            {
                stream_out << t(k, q, alpha) << " ";
            }
            stream_out << std::endl;
        }
    }
}
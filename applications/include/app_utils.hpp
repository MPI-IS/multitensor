// Copyright (c) 2019, Max Planck Society / Software Workshop - Max Planck Institute for Intelligent Systems
// Distributed under the GNU GPL license version 3
// See file LICENSE.md or at https://github.com/MPI-IS/multitensor/LICENSE.md

/*!
 * @file
 *
 * @brief Application utilities.
 *
 * @author Jean-Claude Passy (jean-claude.passy@tuebingen.mpg.de)
 * @author Caterina De Bacco (caterina.debacco@tuebingen.mpg.de)
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

#include "multitensor/params.hpp"
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
 * @tparam weight_t Affinity values type
 *
 * @param[in] filename Name of the file containing the data
 * @param[in,out] edges_start Labels of vertices where an edge starts
 * @param[in,out] edges_end Labels of vertices where an edge ends
 * @param[in,out] edges_weight Edge weights
 */
template <class weight_t>
void read_adjacency_data(const boost::filesystem::path &filename,
                         std::vector<size_t> &edges_start,
                         std::vector<size_t> &edges_end,
                         std::vector<weight_t> &edges_weight)
{
    assert(edges_start.size() == 0);
    assert(edges_end.size() == 0);
    assert(edges_weight.size() == 0);
    std::ifstream in(filename.string());
    if (in.fail())
    {
        throw std::runtime_error(
            std::string("In read_adjacency_data, failed to open ") + filename.string());
    }

    std::cout << "Reading adjacency file " << filename << std::endl;

    std::string line;
    while (!in.eof())
    {
        std::getline(in, line);
        // skip over empty lines
        if (line.size() == 0)
        {
            continue;
        }

        // Remove trailing whitespaces
        line.erase(line.find_last_not_of(" ") + 1);

        std::vector<weight_t> current_weights;
        std::istringstream is(line);

        // First character should be an E for edge
        // std::string tok;
        // is >> tok;
        // assert(tok == "E");

        // Second two characters are the edges ids
        // TO DO: might use tuple later on
        size_t current_edge_in, current_edge_out;
        is >> current_edge_in >> current_edge_out;

        // Read the rest of the data
        weight_t value;
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
 * @brief Read affinity file matrix data and creates the necessary tensor
 *
 * @param[in] filename Name of the file containing the data
 * @param[in] assortative Whether the affinity matrix is assortative
 * @param[in,out] w Vector containing the values of the Affinity tensor
 */
void read_affinity_data(const boost::filesystem::path &filename,
                        const bool &assortative,
                        std::vector<double> &w);

/*!
 * @brief Write affinity file
 *
 * @tparam weight_t Affinity values type
 *
 * @param[in] output_filename Name of the output file
 * @param[in] affinity Vector containing the values of the Affinity tensor
 * @param[in] results Results of the run
 */
template <class weight_t>
void write_affinity_file(const boost::filesystem::path &output_filename,
                         const std::vector<weight_t> &affinity,
                         const utils::Report &results,
                         const size_t nof_groups,
                         const size_t nof_layers)
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
    // Check the dimensions to know if the affinity tensor was assortative
    const size_t assortative = (affinity.size() == (nof_layers * nof_groups));
    size_t index(0);
    for (dimension_t alpha = 0; alpha < nof_layers; alpha++)
    {
        stream_out << "a= " << alpha << std::endl;
        for (dimension_t k = 0; k < nof_groups; k++)
        {
            // Assortative case - we check
            if (assortative)
            {
                index = k + alpha * nof_groups;
                stream_out << affinity[index] << " ";
            }
            else
            {
                for (dimension_t q = 0; q < nof_groups; q++)
                {
                    // The data is transposed (rows are enumerated first)
                    // but we want to go through the columns for a given row
                    // This is why we switch the q and k indices
                    index = k + q * nof_groups + alpha * nof_groups * nof_groups;
                    stream_out << affinity[index] << " ";
                }
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
 * @brief Write a membership matrix into a file
 *
 * @tparam scalar_t Matrix values type
 *
 * @param[in] output_filename Name of the output file
 * @param[in] labels Vertices labels
 * @param[in] mat Matrix (u or v)
 * @param[in] results Results of the run
 */
template <class scalar_t>
void write_membership_file(const boost::filesystem::path &output_filename,
                           const std::vector<size_t> &labels,
                           const tensor::Matrix<scalar_t> &mat,
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
    auto dims = mat.dims();
    const dimension_t nrows{std::get<0>(dims)};
    const dimension_t ncols{std::get<1>(dims)};
    for (size_t k = 0; k < nrows; k++)
    {
        stream_out << labels[k] << " ";
        for (size_t q = 0; q < ncols; q++)
        {
            stream_out << mat(k, q) << " ";
        }
        stream_out << std::endl;
    }
}

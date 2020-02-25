/*!
 * @file
 *
 * @brief Application utilities.
 *
 * @author Jean-Claude Passy (jean-claude.passy@tuebingen.mpg.de)
 */

#pragma once

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include "multitensor_parameters.hpp"

namespace po = boost::program_options;

//! @brief Parsing the command line
std::tuple<po::variables_map, po::options_description> parse_command_line(int ac, char **av);

/*!
 * @brief Read adjacency matrix data and put it in vector
 *
 * @param[in, out] edges Vector containing the edges
 * @param[in, out] edge_weights Vector containing the weight og the edges
 * @param[in, out] nof_nodes Number of nodes
 * @param[in, out] nof_layers Number of layers
 * @param[in] filename Name of the file containing the data
 *
 */

template <class scalar_t>
void read_adjacency_data(std::vector<size_t> &edges_in,
                         std::vector<size_t> &edges_out,
                         std::vector<scalar_t> &edges_weight,
                         size_t &nof_nodes,
                         size_t &nof_layers,
                         const std::string &filename)
{

    edges_in.clear();
    edges_out.clear();
    edges_weight.clear();
    std::ifstream in(filename.c_str());
    if (in.fail())
    {
        throw std::runtime_error(
            std::string("In read_adjacency_data, failed to open ") + filename);
    }

    std::cout << "Reading adjacency file " << filename << std::endl;

    std::string line;
    nof_layers = 0;
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

        // Number of layers should be > 0
        assert(!current_weights.empty());
        if (nof_layers == 0)
        {
            nof_layers = current_weights.size();
        }
        assert(nof_layers == current_weights.size());

        // Insert data at the end of the vectors
        edges_in.push_back(current_edge_in);
        edges_out.push_back(current_edge_out);
        edges_weight.insert(
            std::end(edges_weight), std::begin(current_weights), std::end(current_weights));
    }

    // Close file
    in.close();

    // Calculate number of nodes
    auto min_node_in = std::min_element(edges_in.begin(), edges_in.end());
    auto min_node_out = std::min_element(edges_out.begin(), edges_out.end());
    auto max_node_in = std::max_element(edges_in.begin(), edges_in.end());
    auto max_node_out = std::max_element(edges_out.begin(), edges_out.end());

    nof_nodes = std::max(edges_in[std::distance(edges_in.begin(), max_node_in)],
                         edges_out[std::distance(edges_out.begin(), max_node_out)]) -
                std::min(edges_in[std::distance(edges_in.begin(), min_node_in)],
                         edges_out[std::distance(edges_out.begin(), min_node_out)]) +
                1;
}
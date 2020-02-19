/*!
 * @file
 *
 * @brief Implementation of the multitensor algorithm
 *
 * @author Jean-Claude Passy (jean-claude.passy@tuebingen.mpg.de)
 */

#pragma once

#include <vector>
#include <boost/graph/adjacency_list.hpp>

#include "multitensor_parameters.hpp"
#include "multitensor_graph.hpp"

//! Multitensor namespace
namespace multitensor
{

/*!
 * @brief Multitensor factorization call.
 *
 * Generic algorithm function.
 *
 * @tparam weight_t edges weight type
 *
 */
template <class weight_t>
void multitensor_factorization(const std::vector<size_t> &edges_in,
                               const std::vector<size_t> &edges_out,
                               const std::vector<weight_t> &edges_weight,
                               const size_t &nof_nodes,
                               const size_t &nof_layers)
{
    std::cout << "Inside main call" << std::endl;

    // Check number of nodes
    if (nof_nodes < 2)
    {
        throw std::runtime_error(
            "[multitensor] Number of nodes should be at least 2, intead got " +
            std::to_string(nof_nodes) + "\n");
    }

    // Check number of layers
    if (nof_layers < 1)
    {
        throw std::runtime_error(
            "[multitensor] Number of layers should be at least 1, intead got " +
            std::to_string(nof_layers) + "\n");
    }

    const size_t nof_edges = edges_in.size();

    // Check same size between in and out edges
    if (nof_edges != edges_out.size())
    {
        throw std::runtime_error(
            "[multitensor] Inconsitent edges: " +
            std::to_string(nof_edges) + " source nodes != " +
            std::to_string(edges_out.size()) + " target nodes\n");
    }

    // Check same size with weights
    if (nof_edges != edges_weight.size() / nof_layers)
    {
        throw std::runtime_error(
            "[multitensor] Inconsitent edges: " +
            std::to_string(nof_edges) + " source nodes but " +
            std::to_string(edges_weight.size() / nof_layers) + " weights provided\n");
    }

    // Build network, i.e. a vector of graphs, one for each layer
    std::vector<graph::Graph<>> A(nof_layers), Aold(nof_layers);
    graph::build_network(A, edges_in, edges_out, edges_weight, nof_nodes);

    // Print some information
    graph::print_network_stat(A);
}

} // namespace multitensor
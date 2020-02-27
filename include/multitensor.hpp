/*!
 * @file
 *
 * @brief Implementation of the multitensor algorithm
 *
 * @author Jean-Claude Passy (jean-claude.passy@tuebingen.mpg.de)
 */

#pragma once

#include <vector>
#include <random>
#include <chrono>
#include <boost/graph/adjacency_list.hpp>
#include <boost/random.hpp>

#include "multitensor_parameters.hpp"
#include "multitensor_graph.hpp"
#include "multitensor_solver.hpp"
#include "multitensor_tensor.hpp"

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
template <
    class graph_t = graph::Graph<>,
    class vertex_t,
    class weight_t>
void multitensor_algo(const std::vector<vertex_t> &edges_in,
                      const std::vector<vertex_t> &edges_out,
                      const std::vector<weight_t> &edges_weight,
                      const size_t &nof_nodes,
                      const size_t &nof_layers,
                      const size_t &nof_groups,
                      const size_t &nof_realizations,
                      const size_t &max_nof_iterations,
                      const size_t &nof_convergences)
{
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
    if (nof_edges * nof_layers != edges_weight.size())
    {
        throw std::runtime_error(
            "[multitensor] Inconsitent edges: " +
            std::to_string(nof_edges) + " source nodes but " +
            std::to_string(edges_weight.size() / nof_layers) + " weights provided\n");
    }

    // Check number of groups
    if (nof_groups < 2)
    {
        throw std::runtime_error(
            "[multitensor] Number of groups should be at least 2, intead got " +
            std::to_string(nof_groups) + "\n");
    }

    // Check number of realizations
    if (nof_realizations < 1)
    {
        throw std::runtime_error(
            "[multitensor] Number of realizations should be at least 1, intead got " +
            std::to_string(nof_realizations) + "\n");
    }

    // Some initial diagnostics
    std::cout << "Number of vertices: " << nof_nodes << std::endl;
    std::cout << "Number of layers: " << nof_layers << std::endl;
    std::cout << "Number of groups: " << nof_groups << std::endl;

    // Start timer
    using clock_t = std::chrono::high_resolution_clock;
    clock_t::time_point start_clock = clock_t::now();

    // Build multilayer network, i.e. a vector of graphs, one for each layer
    std::vector<graph_t> A(nof_layers);
    graph::build_network(A, edges_in, edges_out, edges_weight);
    graph::print_graph_stats(A);

    // Extract indices of vertices that have at least one out/in-going edges
    std::vector<size_t> index_vertices_with_out_edges;
    std::vector<size_t> index_vertices_with_in_edges;
    graph::extract_vertices_with_edges(index_vertices_with_out_edges, index_vertices_with_in_edges, A);

    // Random generator
    // TO DO: improve, no real seed!
    boost::mt19937 gen;                               // random real and integer number generators
    boost::uniform_real<double> const uni_dist(0, 1); // uniform distribution for real numbers
    // picks a random real number in (0,1)
    boost::variate_generator<boost::mt19937 &, boost::uniform_real<>> random_generator(gen, uni_dist);

    // Build solver and run
    solver::Solver solver(nof_realizations, max_nof_iterations, nof_convergences);
    solver.run(nof_nodes, nof_groups, nof_layers,
               index_vertices_with_out_edges, index_vertices_with_in_edges,
               A, random_generator);

    // Duration
    std::chrono::duration<double, std::milli> time_spam_ms = clock_t::now() - start_clock;

    // Some last outputs
    graph::print_graph_stats(A);
    std::cout << "Total duration: " << time_spam_ms.count() / 1000. << " s" << std::endl;
}

} // namespace multitensor
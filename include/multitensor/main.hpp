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
#include <string>
#include <cstddef>
#include <set>
#include <boost/random.hpp>
#include <boost/filesystem.hpp>

#include "multitensor/parameters.hpp"
#include "multitensor/graph.hpp"
#include "multitensor/initialization.hpp"
#include "multitensor/solver.hpp"
#include "multitensor/tensor.hpp"
#include "multitensor/utils.hpp"

//! Multitensor namespace
namespace multitensor
{

/*!
 * @brief Multitensor factorization algorithm call.
 *
 * @tparam vertex_t Edges label type
 * @tparam weight_t The type of an edge weight
 *
 * @param[in] edges_start Labels of vertices where an edge starts
 * @param[in] edges_end Labels of vertices where an edge ends
 * @param[in] edges_weight Edges weights
 * @param[in] output_directory Name of the output directory
 * @param[in] nof_groups Number of groups
 * @param[in] nof_realizations Number of realizations
 * @param[in] max_nof_iterations Maximum number of iterations in each realization
 * @param[in] nof_convergences Number of successive passed convergence criteria for declaring the results converged
 *
 * @returns Detailed results of the algorithm.
 */
template <class vertex_t,
          class weight_t,
          class w_init_t = initialization::init_tensor_symmetric_random,
          class uv_init_t = initialization::init_tensor_partial_random>
utils::Report multitensor_factorization(const std::vector<vertex_t> &edges_start,
                                        const std::vector<vertex_t> &edges_end,
                                        const std::vector<weight_t> &edges_weight,
                                        const std::string &output_directory,
                                        const size_t &nof_groups,
                                        const size_t &nof_realizations,
                                        const size_t &max_nof_iterations,
                                        const size_t &nof_convergences)
{
    // Definitions and checks
    // Number of edges
    const size_t nof_edges = edges_start.size();
    if (nof_edges < 1)
    {
        throw std::runtime_error(
            "[multitensor] Number of edges should be at least 1, intead got " +
            std::to_string(nof_edges) + "\n");
    }
    // Check same size between in and out edges
    if (nof_edges != edges_end.size())
    {
        throw std::runtime_error(
            "[multitensor] Inconsitent edges: " +
            std::to_string(nof_edges) + " source vertices != " +
            std::to_string(edges_end.size()) + " target vertices\n");
    }

    // Number of vertices
    const size_t nof_vertices = utils::calculate_num_vertices(edges_start, edges_end);
    if (nof_vertices < 2)
    {
        throw std::runtime_error(
            "[multitensor] Number of vertices should be at least 2, intead got " +
            std::to_string(nof_vertices) + "\n");
    }

    // Number of layers
    const size_t nof_layers = edges_weight.size() / nof_edges;
    if (nof_layers < 1)
    {
        throw std::runtime_error(
            "[multitensor] Number of layers should be at least 1, intead got " +
            std::to_string(nof_layers) + "\n");
    }

    // Number of groups
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

    // Check max number of iterations
    if (max_nof_iterations < 1)
    {
        throw std::runtime_error(
            "[multitensor] Maximum number of iterations should be at least 1, intead got " +
            std::to_string(max_nof_iterations) + "\n");
    }

    // Check number of convergences
    if (nof_convergences < 1)
    {
        throw std::runtime_error(
            "[multitensor] Number of convergences should be at least 1, intead got " +
            std::to_string(nof_convergences) + "\n");
    }

    // Start timer
    using clock_t = std::chrono::high_resolution_clock;
    clock_t::time_point start_clock = clock_t::now();

    // Build multilayer network
    graph::Network A(edges_start, edges_end, edges_weight);
    A.print_graph_stats();

    // Extract indices of vertices that have at least one out/in-going edges
    std::vector<size_t> index_vertices_with_out_edges;
    std::vector<size_t> index_vertices_with_in_edges;
    A.extract_vertices_with_edges(index_vertices_with_out_edges, index_vertices_with_in_edges);
    assert(index_vertices_with_out_edges.size() <= nof_vertices);
    assert(index_vertices_with_in_edges.size() <= nof_vertices);

    // Random generator
    // TO DO: improve, no real seed!
    boost::mt19937 gen;                               // random real and integer number generators
    boost::uniform_real<double> const uni_dist(0, 1); // uniform distribution for real numbers
    // picks a random real number in (0,1)
    boost::variate_generator<boost::mt19937 &, boost::uniform_real<>> random_generator(gen, uni_dist);

    // Affinity tensor for the groups
    tensor::Tensor<double> w(nof_groups, nof_groups, nof_layers);
    // Tensors linking vertices in groups
    tensor::Tensor<double> u(nof_vertices, nof_groups), v(nof_vertices, nof_groups);

    // Run solver
    // Diagnostics
    std::cout << "Starting algorithm...\n"
              << "\t- Number of vertices: " << nof_vertices << " ("
              << index_vertices_with_out_edges.size() << " with outgoing edges, "
              << index_vertices_with_in_edges.size() << " with incoming edges)\n"
              << "\t- Number of layers: " << nof_layers << "\n"
              << "\t- Number of groups: " << nof_groups << "\n"
              << "\t- Number of realizations: " << nof_realizations << "\n"
              << std::endl;
    solver::Solver solver(nof_realizations, max_nof_iterations, nof_convergences);
    utils::Report results = solver.run<w_init_t, uv_init_t>(
        index_vertices_with_out_edges, index_vertices_with_in_edges, A, w, u, v, random_generator);

    // Update duration
    std::chrono::duration<double, std::milli> time_spam_ms = clock_t::now() - start_clock;
    results.duration = time_spam_ms.count() / 1000.;

    // Outputs
    std::cout << "Done!" << std::endl;
    printf("Likelihood: %8.3e\n", results.max_L2());
    printf("Duration: %8.3e s\n", results.duration);

    // Write files - might be moved somewhere else later on
    boost::filesystem::create_directory(output_directory);
    auto dpath = boost::filesystem::canonical(output_directory);
    std::cout << "Writing output files in directory " << dpath << std::endl;
    utils::write_affinity_file(dpath / WOUT_FILENAME, w, results);
    utils::write_membership_file(dpath / UOUT_FILENAME, A, u, results);
    utils::write_membership_file(dpath / VOUT_FILENAME, A, v, results);
    utils::write_info_file(dpath / INFO_FILENAME, results);

    // Returns
    return results;
}
} // namespace multitensor
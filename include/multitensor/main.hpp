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
#include <boost/graph/adjacency_list.hpp>

#include "multitensor/params.hpp"
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
 * @tparam direction_t The type of graph (directed or undirected)
 * @tparam affinity_t The type of tensor for the affinity
 * @tparam affinity_init_t The type of the initialization for the affinity
 * @tparam vertex_t Edges label type
 * @tparam weight_t The type of an edge weight
 * @tparam random_t The type of the random generator
 *
 * @param[in] edges_start Labels of vertices where an edge starts
 * @param[in] edges_end Labels of vertices where an edge ends
 * @param[in] edges_weight Edges weights
 * @param[in] nof_realizations Number of realizations
 * @param[in] max_nof_iterations Maximum number of iterations in each realization
 * @param[in] nof_convergences Number of successive passed convergence criteria for declaring the results converged
 * @param[in,out] labels Vertices labels
 * @param[in,out] u Tensors linking outgoing vertices
 * @param[in,out] v Tensors linking incoming vertices
 * @param[in,out] affinity Affinity values
 * @param[in,out] random_generator Random generator
 *
 * @returns Detailed results of the algorithm.
 */
template <class direction_t = boost::bidirectionalS,
          class affinity_t = tensor::SymmetricTensor<double>,
          class affinity_init_t = initialization::init_symmetric_tensor_random,
          class vertex_t,
          class weight_t,
          class random_t = utils::RandomGenerator<>>
utils::Report multitensor_factorization(const std::vector<vertex_t> &edges_start,
                                        const std::vector<vertex_t> &edges_end,
                                        const std::vector<weight_t> &edges_weight,
                                        const size_t &nof_realizations,
                                        const size_t &max_nof_iterations,
                                        const size_t &nof_convergences,
                                        std::vector<vertex_t> &labels,
                                        tensor::Matrix<double> &u,
                                        tensor::Matrix<double> &v,
                                        std::vector<double> &affinity,
                                        random_t random_generator = random_t{})
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

    // Number of layers
    if (edges_weight.size() % nof_edges != 0)
    {
        throw std::runtime_error(
            "[multitensor] Number of weights should be a multiple of the number of edges, intead got " +
            std::to_string(edges_weight.size()) + " weights and " +
            std::to_string(nof_edges) + " edges\n");
    }
    const size_t nof_layers = edges_weight.size() / nof_edges;
    if (nof_layers < 1)
    {
        throw std::runtime_error(
            "[multitensor] Number of layers should be at least 1, intead got " +
            std::to_string(nof_layers) + "\n");
    }

    // Number of groups
    size_t nof_groups(0), affinity_size(0);
    // Assortative case
    if constexpr (std::is_same_v<affinity_t, tensor::DiagonalTensor<double>>)
    {
        nof_groups = static_cast<size_t>(affinity.size() / nof_layers);
        affinity_size = nof_groups * nof_layers;
    }
    else // Non-assortative case
    {
        nof_groups = static_cast<size_t>(std::sqrt(affinity.size() / nof_layers));
        affinity_size = nof_groups * nof_groups * nof_layers;
    }
    if (nof_groups < 2)
    {
        throw std::runtime_error(
            "[multitensor] Number of groups should be at least 2, intead got " +
            std::to_string(nof_groups) + "\n");
    }
    if (affinity_size != affinity.size())
    {
        throw std::runtime_error(
            "[multitensor] W size should have the form (k x k x nof_layers) with k = " +
            std::to_string(nof_groups) +
            " and nof_layers = " + std::to_string(nof_layers) +
            ", intead got " + std::to_string(affinity.size()) + "\n");
    }

    // Number of vertices
    const size_t nof_vertices = utils::get_num_vertices(edges_start, edges_end);
    if (nof_vertices < 2)
    {
        throw std::runtime_error(
            "[multitensor] Number of vertices should be at least 2, intead got " +
            std::to_string(nof_vertices) + "\n");
    }
    if ((nof_vertices * nof_groups) != u.size())
    {
        throw std::runtime_error(
            "[multitensor] U size should have the form (k x nof_vertices) with k = " +
            std::to_string(nof_groups) +
            " and nof_vertices = " + std::to_string(nof_vertices) +
            ", intead got " + std::to_string(u.size()) + "\n");
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
    graph::Network<vertex_t, direction_t> A(edges_start, edges_end, edges_weight);
    A.print_graph_stats();

    // Extract indices of vertices that have at least one out/in-going edges
    auto u_list = std::make_shared<std::vector<size_t>>();
    auto v_list = std::make_shared<std::vector<size_t>>();
    A.extract_vertices_with_edges(u_list, v_list);

    // Create affinity tensor by copy
    // @todo Later on, we might use a reference to the orginal vector instead
    affinity_t w(nof_groups, nof_layers, affinity);

    // Run solver
    // Diagnostics
    std::cout << "Starting algorithm...\n"
              << "\t- Number of vertices: " << nof_vertices << "\n"
              << "\t- Number of layers: " << nof_layers << "\n"
              << "\t- Number of groups: " << nof_groups << "\n"
              << "\t- Number of realizations: " << nof_realizations << "\n"
              << std::endl;
    solver::Solver solver(nof_realizations, max_nof_iterations, nof_convergences);
    utils::Report results = solver.run<affinity_init_t>(
        *u_list, *v_list, A, u, v, w, random_generator);

    // Copy Network labels for outputs and copy affinity data
    // @todo improve, too much data copy
    A.extract_vertices_labels(labels);
    affinity = w.get_data();

    // Update duration and seed
    std::chrono::duration<double, std::milli> time_spam_ms = clock_t::now() - start_clock;
    results.duration = time_spam_ms.count() / 1000.;
    results.seed = random_generator.seed;

    // Outputs
    std::cout << "Done!" << std::endl;
    printf("Likelihood: %8.3e\n", results.max_L2());
    printf("Duration: %8.3e s\n", results.duration);

    // Returns
    return results;
}
} // namespace multitensor
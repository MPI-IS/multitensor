/*!
 * @file
 *
 * @brief Implementation of graphs and related functions
 *
 * @author Jean-Claude Passy (jean-claude.passy@tuebingen.mpg.de)
 *
 * This file contains the implementation for creating and manipulating
 * graphs, and building a multilayer network of graphs.
 */

#pragma once

#include <iostream>
#include <vector>
#include <boost/graph/adjacency_list.hpp>

namespace multitensor
{

//! Graph and network manipulation
namespace graph
{

/*!
 * @brief Class for the Vertex properties.
 *
 * It allows to distinguish between node id and node label.
 */
struct VertexProperty
{
    //! Node label
    const std::string label;
};

template <class list_t = boost::vecS>
using Graph = boost::adjacency_list<list_t, list_t, boost::bidirectionalS, VertexProperty>;

template <class graph_t = Graph<>>
using Vertex = typename boost::graph_traits<graph_t>::vertex_descriptor;

template <class graph_t = Graph<>>
using out_edge_iterator = typename boost::graph_traits<graph_t>::out_edge_iterator;

template <class graph_t = Graph<>>
using in_edge_iterator = typename boost::graph_traits<graph_t>::in_edge_iterator;

/*!
 * @brief Adding a vertex to a network
 *
 * @tparam graph_t Graph type
 *
 * @param[in] label Node label
 * @param[in,out] A Network
 *
 * @return The node id
 */
template <class graph_t = Graph<>>
size_t add_vertex(const std::string &label, std::vector<graph_t> &A)
{
    static std::map<std::string, size_t> idx_map;
    std::map<std::string, size_t>::iterator it = idx_map.find(label);
    size_t idx;

    if (it == idx_map.end())
    {
        for (size_t alpha = 0; alpha < A.size(); alpha++)
        {

            idx = boost::add_vertex(VertexProperty{label}, A[alpha]);
        }
        idx_map[label] = idx;
        return idx_map[label];
    }
    return it->second;
}

/*!
 * @brief Printing information about the network
 *
 * @tparam graph_t Graph type
 *
 * @param[in] A Network
 */
template <class graph_t>
void print_graph_stats(std::vector<graph_t> &A)
{
    size_t nof_layers = A.size();
    size_t num_edges_in_layer; // number of edges in each layer
    size_t N = boost::num_vertices(A[0]);
    double density;

    for (size_t alpha = 0; alpha < nof_layers; alpha++)
    {
        num_edges_in_layer = boost::num_edges(A[alpha]);
        density = 100. * (double)num_edges_in_layer / (double)(N * (N - 1.));
        std::cout << "E[" << alpha << "] = " << num_edges_in_layer
                  << "  density= " << density
                  << std::endl;
    }
}

/*!
 * @brief Building a multilayer network.
 *
 * @tparam graph_t Graph type
 * @tparam vertex_t Edges label type
 * @tparam weight_t The type of an edge weight
 *
 * @param[in] edges_in Labels of in-going vertices for each edge
 * @param[in] edges_out Labels of out-going vertices for each edge
 * @param[in] edges_weight Edge weights
 * @param[in,out] A Network
 */
template <class graph_t,
          class vertex_t,
          class weight_t>
void build_network(const std::vector<vertex_t> &edges_in,
                   const std::vector<vertex_t> &edges_out,
                   const std::vector<weight_t> &edges_weight,
                   std::vector<graph_t> &A)
{
    const size_t nof_layers = A.size();

    for (size_t i = 0; i < edges_in.size(); i++)
    {
        // Adding vertices
        size_t node_in = add_vertex(std::to_string(edges_in[i]), A);
        size_t node_out = add_vertex(std::to_string(edges_out[i]), A);

        // Adding edges if needed
        for (size_t alpha = 0; alpha < nof_layers; alpha++)
        {
            weight_t weight = edges_weight[i * nof_layers + alpha];
            if (weight > EPS_PRECISION)
            {
                for (size_t w = 0; w < weight; w++) // seems dangerous...
                {
                    boost::add_edge(
                        boost::vertex(node_in, A[alpha]),
                        boost::vertex(node_out, A[alpha]),
                        A[alpha]);
                }
            }
        }
    }
}

/*!
 * @brief Extracting vertices with out/in going edges.
 *
 * @tparam graph_t graph type
 *
 * @param[in] A Netowrk
 * @param[in, out] index_vertices_with_out_edges Indices of the vertices with at least one outgoing edege
 * @param[in, out] index_vertices_with_in_edges Indices of the vertices with at least one incoming edege
 */
template <class graph_t>
void extract_vertices_with_edges(const std::vector<graph_t> &A,
                                 std::vector<size_t> &index_vertices_with_out_edges,
                                 std::vector<size_t> &index_vertices_with_in_edges)
{
    unsigned int nof_edges_out, nof_edges_in;
    out_edge_iterator<graph_t> eit_out, eend_out;
    in_edge_iterator<graph_t> eit_in, eend_in;

    // Loop through all the vertices in the first layer
    for (size_t i = 0; i < boost::num_vertices(A[0]); i++)
    {
        nof_edges_out = nof_edges_in = 0;

        for (size_t alpha = 0; alpha < A.size(); alpha++)
        {
            // Number of edges out
            for (std::tie(eit_out, eend_out) = boost::out_edges(i, A[alpha]); eit_out != eend_out; ++eit_out)
            {
                nof_edges_out++;
            }

            // Number of edges in
            for (std::tie(eit_in, eend_in) = boost::in_edges(i, A[alpha]); eit_in != eend_in; ++eit_in)
            {
                nof_edges_in++;
            }
        }

        // Record edges if necessary
        if (nof_edges_out)
        {
            index_vertices_with_out_edges.emplace_back(i);
        }
        if (nof_edges_in)
        {
            index_vertices_with_in_edges.emplace_back(i);
        }
    }

    // Diagnostics
    std::cout << "Number of vertices with outgoing edges: " << index_vertices_with_out_edges.size() << std::endl;
    std::cout << "Number of vertices with ingoing edges: " << index_vertices_with_in_edges.size() << std::endl;
}

} // namespace graph
} // namespace multitensor
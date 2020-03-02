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
    std::string label;
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
 * @brief Class representing a multilayer network
 */
template <class graph_t = Graph<>>
class Network
{
private:
    dimension_t nlayers;
    unsigned int nvertices, nedges;
    std::vector<graph_t> layers;

    /*!
     * @brief Adding a vertex to a network
     */
    template <class vertex_t>
    size_t add_vertex(const vertex_t &label)
    {
        static std::map<vertex_t, size_t> idx_map;
        typename std::map<vertex_t, size_t>::iterator it = idx_map.find(label);
        size_t idx(0);

        if (it == idx_map.end())
        {
            for (size_t alpha = 0; alpha < nlayers; alpha++)
            {
                idx = boost::add_vertex(VertexProperty{label}, operator()(alpha));
            }
            idx_map[label] = idx;
            return idx_map[label];
        }
        return it->second;
    }

public:
    /*!
     * @brief Network constructor
     */
    template <class vertex_t,
              class weight_t>
    Network(const std::vector<vertex_t> &edges_in,
            const std::vector<vertex_t> &edges_out,
            const std::vector<weight_t> &edges_weight)
        : nlayers((edges_in.size() == 0) ? 0 : edges_weight.size() / edges_in.size()),
          nvertices(0),
          nedges(0),
          layers(nlayers, graph_t(0))
    {

        assert(edges_in.size() == edges_out.size());
        assert((nlayers * edges_in.size()) == edges_weight.size());
        for (size_t i = 0; i < edges_in.size(); i++)
        {
            // Adding vertices
            size_t node_in = add_vertex(std::to_string(edges_in[i]));
            size_t node_out = add_vertex(std::to_string(edges_out[i]));

            // Adding edges if needed
            for (size_t alpha = 0; alpha < nlayers; alpha++)
            {
                weight_t weight = edges_weight[i * nlayers + alpha];
                if (weight > EPS_PRECISION)
                {
                    for (size_t w = 0; w < weight; w++) // seems dangerous...
                    {
                        boost::add_edge(boost::vertex(node_in, operator()(alpha)),
                                        boost::vertex(node_out, operator()(alpha)),
                                        operator()(alpha));
                        nedges++;
                    }
                }
            }
        }
        // Store number of vertices
        if (nedges > 0)
        {
            nvertices = boost::num_vertices(operator()(0));
        }
    }

    /*!
     * @brief Extracting vertices with out/in going edges.
     *
     * @param[in, out] index_vertices_with_out_edges Indices of the vertices with at least one outgoing edege
     * @param[in, out] index_vertices_with_in_edges Indices of the vertices with at least one incoming edege
     */
    void extract_vertices_with_edges(std::vector<size_t> &index_vertices_with_out_edges,
                                     std::vector<size_t> &index_vertices_with_in_edges)
    {
        unsigned int nof_edges_out, nof_edges_in;
        out_edge_iterator<graph_t> eit_out, eend_out;
        in_edge_iterator<graph_t> eit_in, eend_in;

        // Loop through all the vertices in the first layer
        for (size_t i = 0; i < num_vertices(); i++)
        {
            nof_edges_out = nof_edges_in = 0;

            for (size_t alpha = 0; alpha < size(); alpha++)
            {
                // Number of edges out
                for (std::tie(eit_out, eend_out) = boost::out_edges(i, operator()(alpha)); eit_out != eend_out; ++eit_out)
                {
                    nof_edges_out++;
                }

                // Number of edges in
                for (std::tie(eit_in, eend_in) = boost::in_edges(i, operator()(alpha)); eit_in != eend_in; ++eit_in)
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

    /*!
     * @brief Returns a reference to a layer in the network
     */
    graph_t &operator()(const size_t alpha)
    {
        return layers[alpha];
    }

    /*!
     * @brief Returns a const reference to a layer in the network
     */
    const graph_t &operator()(const size_t alpha) const
    {
        return layers[alpha];
    }

    //! @brief Returns the number of layers in the network
    size_t size() const noexcept
    {
        return nlayers;
    }

    //! @brief Returns the number of vertices in the network
    size_t num_vertices() const noexcept
    {
        return nvertices;
    }

    //! @brief Returns the number of edges in the full network
    size_t num_edges() const noexcept
    {
        return nedges;
    }

    //! @brief Printing information about the network
    void print_graph_stats()
    {
        size_t num_edges_in_layer;
        size_t N = num_vertices();
        double density;
        std::cout << "N = " << N << std::endl;

        for (size_t alpha = 0; alpha < size(); alpha++)
        {
            num_edges_in_layer = boost::num_edges(operator()(alpha));
            density = 100. * (double)num_edges_in_layer / (double)(N * (N - 1.));
            std::cout << "E[" << alpha << "] = " << num_edges_in_layer
                      << "  density= " << density
                      << std::endl;
        }
    }
};

} // namespace graph
} // namespace multitensor
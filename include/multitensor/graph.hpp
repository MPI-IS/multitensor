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
#include <cstddef>
#include <boost/graph/adjacency_list.hpp>

namespace multitensor
{

//! Graph and network manipulation
namespace graph
{

/*!
 * @brief Class for the Vertex properties.
 *
 * @tparam vertex_t Vertex label type
 *
 */
template <class vertex_t = std::string>
struct VertexProperty
{
    //! Node label
    vertex_t label;
};

template <class list_t = boost::vecS>
using Graph = boost::adjacency_list<list_t, list_t, boost::bidirectionalS, VertexProperty<>>;

template <class graph_t = Graph<>>
using Vertex = typename boost::graph_traits<graph_t>::vertex_descriptor;

template <class graph_t = Graph<>>
using out_edge_iterator = typename boost::graph_traits<graph_t>::out_edge_iterator;

template <class graph_t = Graph<>>
using in_edge_iterator = typename boost::graph_traits<graph_t>::in_edge_iterator;

/*!
 * @brief Class representing a multilayer network
 *
 * @tparam direction_t Graph type (directed or undirected)
 * @tparam vertex_t Vertex label type
 */
template <class vertex_t,
          class direction_t = boost::bidirectionalS>
class Network
{
    using graph_t = boost::adjacency_list<boost::vecS, boost::vecS, direction_t, VertexProperty<vertex_t>>;
    using out_edge_iterator_t = typename boost::graph_traits<graph_t>::out_edge_iterator;
    using in_edge_iterator_t = typename boost::graph_traits<graph_t>::in_edge_iterator;

private:
    dimension_t nlayers;
    unsigned int nvertices, nedges;
    std::vector<graph_t> layers;

    //! @brief Mapping vertex label <-> index
    std::map<vertex_t, size_t> idx_map;

    /*!
     * @brief Adding a vertex to a network
     *
     * @param[in] label Vertex label
     *
     * @returns Vertex index
     *
     * @note Vertices are not duplicated
     */
    size_t add_vertex(const vertex_t &label)
    {
        typename std::map<vertex_t, size_t>::iterator it = idx_map.find(label);
        size_t idx(0);

        if (it == idx_map.end())
        {
            for (size_t alpha = 0; alpha < nlayers; alpha++)
            {
                idx = boost::add_vertex(VertexProperty<vertex_t>{label}, operator()(alpha));
            }
            idx_map[label] = idx;
            return idx_map[label];
        }
        return it->second;
    }

public:
    using direction_type = direction_t;

    /*!
     * @brief Network constructor
     *
     * param[in] edges_start Labels of vertices where an edge starts
     * param[in] edges_end Labels of vertices where an edge ends
     * param[in] edges_weight Edges weights
     *
     * @pre @c edges_start.size() == edges_end.size()
     * @pre @c edges_weight.size() == num_layers * edges_end.size()
     */
    template <class weight_t>
    Network(const std::vector<vertex_t> &edges_start,
            const std::vector<vertex_t> &edges_end,
            const std::vector<weight_t> &edges_weight)
        : nlayers((edges_start.size() == 0) ? 0 : edges_weight.size() / edges_start.size()),
          nvertices(0),
          nedges(0),
          layers(nlayers, graph_t(0))
    {

        assert(edges_start.size() == edges_end.size());
        assert((nlayers * edges_start.size()) == edges_weight.size());
        for (size_t i = 0; i < edges_start.size(); i++)
        {
            // Adding vertices
            size_t node_in = add_vertex(edges_start[i]);
            size_t node_out = add_vertex(edges_end[i]);

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
        if (nlayers > 0)
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
    void extract_vertices_with_edges(std::shared_ptr<std::vector<size_t>> &index_vertices_with_out_edges,
                                     std::shared_ptr<std::vector<size_t>> &index_vertices_with_in_edges)
    {
        // Iterators
        std::pair<out_edge_iterator_t, out_edge_iterator_t> its_out;
        std::pair<in_edge_iterator_t, in_edge_iterator_t> its_in;

        // Loop through all the vertices in the first layer
        for (size_t i = 0; i < num_vertices(); i++)
        {
            // Outgoing edges
            for (size_t alpha = 0; alpha < num_layers(); alpha++)
            {
                its_out = boost::out_edges(i, operator()(alpha));
                if (std::get<0>(its_out) != std::get<1>(its_out))
                {
                    (*index_vertices_with_out_edges).emplace_back(i);
                    break;
                }
            }

            // Compile-time if statement
            // We count incoming edges only if we are using a directed network
            if constexpr (std::is_same_v<direction_t, boost::bidirectionalS>)
            {
                // Incoming edges
                for (size_t alpha = 0; alpha < num_layers(); alpha++)
                {
                    its_in = boost::in_edges(i, operator()(alpha));
                    if (std::get<0>(its_in) != std::get<1>(its_in))
                    {
                        (*index_vertices_with_in_edges).emplace_back(i);
                        break;
                    }
                }
            } // end if consexpr
        }     // end loop vertices

        // Reset pointer if graphs are undirected
        // and check number of shared instances
        if constexpr (std::is_same_v<direction_t, boost::bidirectionalS>)
        {
            assert(!index_vertices_with_out_edges ||
                   index_vertices_with_out_edges.use_count() == 1);
            assert(!index_vertices_with_out_edges ||
                   index_vertices_with_in_edges.use_count() == 1);
        }
        else
        {
            index_vertices_with_in_edges.reset();
            index_vertices_with_in_edges = index_vertices_with_out_edges;
            assert(!index_vertices_with_out_edges ||
                   index_vertices_with_out_edges.use_count() == 2);
            assert(!index_vertices_with_out_edges ||
                   index_vertices_with_in_edges.use_count() == 2);
        }
    }

    /*!
     * @brief Extracting vertices labels.
     *
     * @param[in, out] labels Vertices labels
     */
    void extract_vertices_labels(std::vector<size_t> &labels)
    {
        labels.clear();
        for (size_t i = 0; i < num_vertices(); i++)
        {
            labels.emplace_back(operator()(0)[i].label);
        }
    }

    //! @brief Returns a reference to a layer in the network
    graph_t &operator()(const size_t alpha)
    {
        return layers[alpha];
    }

    //! @brief Returns a const reference to a layer in the network
    const graph_t &operator()(const size_t alpha) const
    {
        return layers[alpha];
    }

    //! @brief Returns the number of layers in the network
    size_t num_layers() const noexcept
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
        size_t num_edges_start_layer;
        size_t N = num_vertices();
        double density;
        std::cout << "N = " << N << std::endl;

        for (size_t alpha = 0; alpha < num_layers(); alpha++)
        {
            num_edges_start_layer = boost::num_edges(operator()(alpha));
            density = 100. * (double)num_edges_start_layer / (double)(N * (N - 1.));
            std::cout << "E[" << alpha << "] = " << num_edges_start_layer
                      << "  density= " << density
                      << std::endl;
        }
    }
};

} // namespace graph
} // namespace multitensor
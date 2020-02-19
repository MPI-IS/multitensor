/*!
 * @file
 *
 * @brief Implementation of a Graph and its capabilities
 *
 * @author Jean-Claude Passy (jean-claude.passy@tuebingen.mpg.de)
 */

#pragma once

#include <vector>
#include <boost/graph/adjacency_list.hpp>

namespace multitensor
{
namespace graph
{

/*!
 * @brief Class for the Vertex properties.
 *
 * It allows to distinguish between node id and node label.
 */
struct VertexProperty
{
    const std::string label; // node label, to be distinguished from the node id
};

template <class list_t = boost::vecS>
using Graph = boost::adjacency_list<list_t, list_t, boost::bidirectionalS, VertexProperty>;

//! @brief Utility function printing some information about the network
template <class graph_t = Graph<>>
void print_network_stat(const std::vector<graph_t> &A)
{
    std::vector<size_t> E(A.size(), 0); // number of edges in each layer
    int N = (int)boost::num_vertices(A[0]);

    for (size_t alpha = 0; alpha < A.size(); alpha++)
    {
        E[alpha] = (size_t)num_edges(A[alpha]);
        std::cout << "E[" << alpha << "] = " << E[alpha]
                  << "  density= " << 100. * ((double)E[alpha] / ((double)(N * (N - 1.))))
                  << std::endl;
    }
}

//! @brief Utility function to add vertex
template <class graph_t = Graph<>>
size_t add_vertex(const std::string &label, std::vector<graph_t> &A)
{
    static std::map<std::string, size_t> idx_map;
    std::map<std::string, size_t>::iterator it = idx_map.find(label);

    if (it == idx_map.end())
    {
        for (size_t alpha = 0; alpha < A.size(); alpha++)
        {
            if (alpha == 0)
            {
                idx_map[label] = boost::add_vertex(VertexProperty{label}, A[alpha]);
            }
            else
            {
                boost::add_vertex(VertexProperty{label}, A[alpha]);
            }
        }
        return idx_map[label];
    }
    return it->second;
}

/*!
 * @brief Building network.
 *
 * @tparam graph_t graph type
 * @tparam weight_t edges weight type
 *
 */
template <class graph_t = Graph<>,
          class weight_t>
void build_network(std::vector<graph_t> &A,
                   const std::vector<size_t> &edges_in,
                   const std::vector<size_t> &edges_out,
                   const std::vector<weight_t> &edges_weight,
                   const size_t &nof_nodes)
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

} // namespace graph
} // namespace multitensor
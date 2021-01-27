// Copyright (c) 2019, Max Planck Society / Software Workshop - Max Planck Institute for Intelligent Systems
// Distributed under the GNU GPL license version 3
// See file LICENSE.md or at https://github.com/MPI-IS/multitensor/LICENSE.md

/*!
 * @file
 *
 * @brief Testing the Graph class.
 *
 * @author Jean-Claude Passy (jean-claude.passy@tuebingen.mpg.de)
 */

#include <iostream>
#include <random>
#include <ctime>
#include <vector>
#include <set>
#include <cstdlib>
#include <boost/test/unit_test.hpp>
#include <boost/graph/adjacency_list.hpp>

#include "multitensor/params.hpp"
#include "multitensor/graph.hpp"
#include "fixtures.hpp"

using namespace multitensor::graph;

BOOST_FIXTURE_TEST_SUITE(tests_graph, fixture_global)

// Checks the VertexProperty class
BOOST_AUTO_TEST_CASE(test_vertex_properties)
{
    // Strings
    std::string label = "Some random label" + std::to_string(rng());
    VertexProperty<> prop{label};
    BOOST_TEST(prop.label == label);

    // Unsigned ints
    size_t label2 = static_cast<size_t>(rng() % 1000);
    VertexProperty<size_t> prop2{label2};
    BOOST_TEST(prop2.label == label2);

    // Ints
    int label3 = static_cast<int>(rng() % 1000) - 500;
    VertexProperty<int> prop3{label3};
    BOOST_TEST(prop3.label == label3);

    // Doubles
    double label4 = static_cast<double>(rng() % 1000) - 500;
    VertexProperty<double> prop4{label4};
    BOOST_TEST(prop4.label == label4);
}

// Checks implementation for an empty network
BOOST_AUTO_TEST_CASE(test_empty_network)
{
    // Building
    std::vector<size_t> vec_empty{};
    BOOST_TEST(vec_empty.size() == 0);
    Network A(vec_empty, vec_empty, vec_empty);

    BOOST_TEST(A.num_layers() == 0);
    BOOST_TEST(A.num_vertices() == 0);
    BOOST_TEST(A.num_edges() == 0);
    BOOST_CHECK_NO_THROW(A.print_graph_stats());

    // Extracting vertices with edges
    A.extract_vertices_with_edges(u_list, v_list);
    BOOST_TEST((*u_list).size() == 0);
    BOOST_TEST((*v_list).size() == 0);
}

// Checks implementation for a network w/o any edge
BOOST_AUTO_TEST_CASE(test_network_wo_edges)
{
    // Building
    std::vector<size_t> empty_weights(edges_weight.size());
    Network A(edges_start, edges_end, empty_weights);

    BOOST_TEST(A.num_layers() == nof_layers);
    BOOST_TEST(A.num_vertices() == nof_vertices);
    BOOST_TEST(A.num_edges() == 0);
    BOOST_CHECK_NO_THROW(A.print_graph_stats());

    // Extracting vertices with edges
    A.extract_vertices_with_edges(u_list, v_list);
    BOOST_TEST((*u_list).size() == 0);
    BOOST_TEST((*v_list).size() == 0);
}

// Checks implementation for a full network
BOOST_AUTO_TEST_CASE(test_full_network)
{
    // Network
    Network A(edges_start, edges_end, edges_weight);

    BOOST_TEST(A.num_layers() == nof_layers);
    BOOST_TEST(A.num_vertices() == nof_vertices);
    BOOST_TEST(A.num_edges() == nof_edges);
    BOOST_CHECK_NO_THROW(A.print_graph_stats());

    // Check vertices
    // Between layers
    std::vector<std::set<size_t>> A_vertices_labels(A.num_layers());
    for (size_t i = 0; i < boost::num_vertices(A(0)); i++)
    {
        for (size_t alpha = 0; alpha < A.num_layers(); alpha++)
        {
            A_vertices_labels[alpha].insert(A(alpha)[i].label);
        }
    }
    // With expected results
    for (size_t alpha = 0; alpha < A.num_layers(); alpha++)
    {
        BOOST_TEST(A_vertices_labels[alpha] == A_vertices_labels[0]);
    }
    BOOST_TEST(A_vertices_labels[0] == v);

    // Check edges
    // First build the mapping from the theoretical data
    using map_t = std::map<std::tuple<size_t, size_t>, size_t>;
    std::vector<map_t> vec_of_maps(A.num_layers());
    for (size_t i = 0; i < edges_start.size(); i++)
    {
        auto edge = std::make_tuple(edges_start[i], edges_end[i]);
        for (size_t alpha = 0; alpha < nof_layers; alpha++)
        {
            auto w = edges_weight[i * nof_layers + alpha];
            if (w > EPS_PRECISION)
            {
                vec_of_maps[alpha][edge] = w;
            }
        }
    }

    // Now the same from the network + check
    for (size_t alpha = 0; alpha < A.num_layers(); alpha++)
    {
        auto es = boost::edges(A(alpha));
        map_t edges_weight_map;
        map_t::iterator map_it;
        for (auto it = es.first; it != es.second; it++)
        {
            auto in = boost::source(*it, A(alpha));
            auto out = boost::target(*it, A(alpha));
            auto edge = std::make_tuple(A(alpha)[in].label,
                                        A(alpha)[out].label);
            map_it = edges_weight_map.find(edge);
            if (map_it == edges_weight_map.end())
            {
                edges_weight_map[edge] = 1;
            }
            else
            {
                edges_weight_map[edge]++;
            }
        }

        // Mapping should be the same
        BOOST_TEST(edges_weight_map == vec_of_maps[alpha]);
    }

    // Extracting vertices with edges
    A.extract_vertices_with_edges(u_list, v_list);
    BOOST_TEST((*u_list).size() == u_labels_theo.size());
    BOOST_TEST((*v_list).size() == v_labels_theo.size());

    // Compare lists
    // NB: u_labels_theo contains labels, u_list contains indices
    std::vector<size_t> u_labels, v_labels;
    for (auto i : *u_list)
    {
        u_labels.emplace_back(A(0)[i].label);
    }
    std::sort(u_labels.begin(), u_labels.end());
    for (auto i : *v_list)
    {
        v_labels.emplace_back(A(0)[i].label);
    }
    std::sort(v_labels.begin(), v_labels.end());
    BOOST_TEST(u_labels == u_labels_theo);
    BOOST_TEST(v_labels == v_labels_theo);

    // Extract labels
    std::vector<size_t> all_labels, all_labels_theo;
    for (size_t i = 0; i < A.num_vertices(); i++)
    {
        all_labels_theo.emplace_back(A(0)[i].label);
    }
    A.extract_vertices_labels(all_labels);
    BOOST_TEST(all_labels == all_labels_theo);
}

BOOST_AUTO_TEST_SUITE_END()

/*!
 * @file
 *
 * @brief Testing the Graph class.
 *
 * @author Jean-Claude Passy (jean-claude.passy@tuebingen.mpg.de)
 */

#include <boost/test/unit_test.hpp>
#include <iostream>
#include <random>
#include <ctime>
#include <vector>
#include <set>
#include <boost/graph/adjacency_list.hpp>

#include "multitensor/parameters.hpp"
#include "multitensor/graph.hpp"

using namespace multitensor::graph;

// Fixture
struct fixture_graph
{
    std::mt19937_64 rng;
    std::time_t seed;
    size_t nof_layers, nof_vertices, nof_data_entries, nof_edges;
    std::vector<unsigned int> edges_in, edges_out, edges_weight;
    std::set<unsigned int> v;

    fixture_graph()
        : seed(1583326877),
          nof_edges(0)
    {
        BOOST_TEST_MESSAGE("In fixture, the seed is " << seed);
        rng.seed(seed);

        nof_layers = size_t(rng() % 10) + 1;
        nof_data_entries = size_t(rng() % 100) + 2;
        edges_in.reserve(nof_data_entries);
        edges_out.reserve(nof_data_entries);
        edges_weight.reserve(nof_data_entries * nof_layers);
        unsigned int tmp_weight;

        // Fill edges vectors and weights
        for (size_t i = 0; i < nof_data_entries; i++)
        {
            edges_in.emplace_back(static_cast<unsigned int>(rng() % 1000));
            edges_out.emplace_back(static_cast<unsigned int>(rng() % 1000));
            for (size_t alpha = 0; alpha < nof_layers; alpha++)
            {
                tmp_weight = static_cast<unsigned int>(rng() % 100);
                edges_weight.emplace_back(tmp_weight);
                if (tmp_weight > EPS_PRECISION)
                {
                    nof_edges += tmp_weight; // 1 edge per weight unit
                }
            }
        }

        // Count number of vertices
        std::set<unsigned int> v_in(edges_in.begin(), edges_in.end());
        std::set<unsigned int> v_out(edges_out.begin(), edges_out.end());
        std::merge(v_in.begin(), v_in.end(),
                   v_out.begin(), v_out.end(),
                   std::inserter(v, v.begin()));
        nof_vertices = v.size();
    }
};

BOOST_FIXTURE_TEST_SUITE(tests_graph, fixture_graph)

// Checks the VertexProperty class
BOOST_AUTO_TEST_CASE(test_vertex_properties)
{
    std::string label = "Some random label" + std::to_string(rng());
    VertexProperty prop{label};
    BOOST_TEST(prop.label == label);
}

// Checks creating an empty network
BOOST_AUTO_TEST_CASE(test_empty_network)
{
    std::vector<unsigned int> vec_empty{};
    BOOST_TEST(vec_empty.size() == 0);

    Network<> A_empty(vec_empty, vec_empty, vec_empty);
    BOOST_TEST(A_empty.size() == 0);
    BOOST_TEST(A_empty.num_edges() == 0);
    BOOST_TEST(A_empty.num_vertices() == 0);
    BOOST_CHECK_NO_THROW(A_empty.print_graph_stats());
}

// Checks creating a full network
BOOST_AUTO_TEST_CASE(test_init_network)
{
    // Network
    Network<> A(edges_in, edges_out, edges_weight);
    BOOST_TEST(A.size() == nof_layers);
    BOOST_TEST(A.num_vertices() == nof_vertices);
    BOOST_TEST(A.num_edges() == nof_edges);
    BOOST_CHECK_NO_THROW(A.print_graph_stats());

    // Check vertices
    // Between layers
    std::vector<std::set<unsigned int>> A_vertices_labels(A.size());
    for (size_t i = 0; i < boost::num_vertices(A(0)); i++)
    {
        for (size_t alpha = 0; alpha < A.size(); alpha++)
        {
            A_vertices_labels[alpha].insert(std::stoul(A(alpha)[i].label));
        }
    }
    // With expected results
    for (size_t alpha = 0; alpha < A.size(); alpha++)
    {
        BOOST_TEST(A_vertices_labels[alpha] == A_vertices_labels[0]);
    }
    BOOST_TEST(A_vertices_labels[0] == v);

    // Check edges
    // First build the mapping from the theoretical data
    using map_t = std::map<std::tuple<unsigned int, unsigned int>, unsigned int>;
    std::vector<map_t> vec_of_maps(A.size());
    for (size_t i = 0; i < edges_in.size(); i++)
    {
        auto edge = std::make_tuple(edges_in[i], edges_out[i]);
        for (size_t alpha = 0; alpha < nof_layers; alpha++)
        {
            unsigned w = edges_weight[i * nof_layers + alpha];
            if (w > EPS_PRECISION)
            {
                vec_of_maps[alpha][edge] = w;
            }
        }
    }

    // Now the same from the network + check
    for (size_t alpha = 0; alpha < A.size(); alpha++)
    {
        auto es = boost::edges(A(alpha));
        map_t edges_weight_map;
        map_t::iterator map_it;
        for (auto it = es.first; it != es.second; it++)
        {
            auto in = boost::source(*it, A(alpha));
            auto out = boost::target(*it, A(alpha));
            auto edge = std::make_tuple(std::stoul(A(alpha)[in].label),
                                        std::stoul(A(alpha)[out].label));
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
}

BOOST_AUTO_TEST_SUITE_END()

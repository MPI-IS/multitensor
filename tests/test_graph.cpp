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

#include "multitensor_parameters.hpp"
#include "multitensor_graph.hpp"

using namespace multitensor::graph;

// Fixture
struct fixture_graph
{
    std::mt19937_64 rng;
    std::time_t seed;
    size_t nof_layers, nof_vertices, nof_edges;
    std::vector<Graph<>> A;

    size_t nof_data_entries;
    std::vector<unsigned int> edges_in, edges_out, edges_weight;

    fixture_graph()
        : seed(std::time(nullptr)),
          nof_edges(0)
    {
        BOOST_TEST_MESSAGE("In fixture, the seed is " << seed);
        rng.seed(seed);
        nof_layers = size_t(rng() % 10) + 1;
        A = std::vector<Graph<>>(nof_layers);

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
        std::set<unsigned int> v;
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

// Checks the function adding a vertex
BOOST_AUTO_TEST_CASE(test_add_vertex)
{
    std::vector<Graph<>> A(nof_layers);
    BOOST_TEST(A.size() == nof_layers);

    // Label
    std::string label = "id " + std::to_string(size_t(rng() % 10));

    size_t x = add_vertex(label, A);
    for (size_t alpha = 0; alpha < A.size(); alpha++)
    {
        BOOST_TEST(num_vertices(A[alpha]) == 1);
        BOOST_TEST(A[alpha][0].label == label);
    }
}

// Checks the function building the multilayer network
BOOST_AUTO_TEST_CASE(test_build_network)
{

    // Build network
    build_network(A, edges_in, edges_out, edges_weight);

    size_t n_edges(0);
    for (size_t alpha = 0; alpha < nof_layers; alpha++)
    {
        n_edges += boost::num_edges(A[alpha]);
    }
    BOOST_TEST(nof_edges == n_edges);
    BOOST_TEST(nof_vertices == boost::num_vertices(A[0]));
}

// Checks the function extracting vertices with edges
// TO DO: finish this test
BOOST_AUTO_TEST_CASE(test_extract_vertices_with_edges)
{
    // Build network
    build_network(A, edges_in, edges_out, edges_weight);

    // Extract vertices of interest
    std::vector<size_t> u, v;
    extract_vertices_with_edges(u, v, A);

    // Build edges sets to compare
    std::vector<unsigned int> u_compare, v_compare;
    bool edge_exists;

    for (size_t i = 0; i < boost::num_vertices(A[0]); i++)
    {
        edge_exists = false;

        // Find vertex label

        for (size_t alpha = 0; alpha < nof_layers; alpha++)
        {
            if (edges_weight[i * nof_layers + alpha] > EPS_PRECISION)
            {
                edge_exists = true;
                break;
            }
        }

        // Add vertices if edge exists
        if (edge_exists)
        {
            u_compare.push_back(i);
            v_compare.push_back(i);
        }
    }

    // Transform to sets
    std::set<unsigned int> set_u(u_compare.begin(), u_compare.end());
    std::set<unsigned int> set_v(v_compare.begin(), v_compare.end());
    std::set<unsigned int> vertices;

    // Comparisons
    BOOST_TEST(1 == 1);
}

BOOST_AUTO_TEST_SUITE_END()

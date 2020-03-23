/*!
 * @file
 *
 * @brief Test fixtures
 *
 * @author Jean-Claude Passy (jean-claude.passy@tuebingen.mpg.de)
 */

#pragma once

#include <iostream>
#include <random>
#include <vector>
#include <set>
#include <cstdlib>
#include <boost/test/unit_test.hpp>

// Fixture for the random generator
struct fixture_rng
{
    std::mt19937_64 rng;
    std::time_t seed;

    fixture_rng()
        : seed(std::time(nullptr))
    {
        BOOST_TEST_MESSAGE("In fixture_rng, the seed is " << seed);
        rng.seed(seed);
    }
};

// Global fixture preparing data
struct fixture_global : fixture_rng
{
    size_t nof_layers, nof_vertices, nof_data_entries, nof_edges;
    std::vector<unsigned int>
        edges_start, edges_end, edges_weight,
        u_labels_theo, v_labels_theo;
    std::set<unsigned int> v;
    std::shared_ptr<std::vector<size_t>> u_list = std::make_shared<std::vector<size_t>>();
    std::shared_ptr<std::vector<size_t>> v_list = std::make_shared<std::vector<size_t>>();

    fixture_global()
        : nof_layers(rng() % 10 + 1),
          nof_data_entries(rng() % 100 + 2),
          nof_edges(0)
    {
        unsigned int estart, eend, tmp_weight;
        std::vector<unsigned int> u_labels_tmp, v_labels_tmp;

        // Fill edges vectors and weights
        for (size_t i = 0; i < nof_data_entries; i++)
        {
            estart = rng() % 1000;
            eend = rng() % 1000;
            edges_start.emplace_back(estart);
            edges_end.emplace_back(eend);
            for (size_t alpha = 0; alpha < nof_layers; alpha++)
            {
                tmp_weight = rng() % 100;
                edges_weight.emplace_back(tmp_weight);
                if (tmp_weight > EPS_PRECISION)
                {
                    u_labels_tmp.emplace_back(estart);
                    v_labels_tmp.emplace_back(eend);
                    nof_edges += tmp_weight; // 1 edge per weight unit
                }
            }
        }

        // Count number of vertices
        std::set<unsigned int> v_in(edges_start.begin(), edges_start.end());
        std::set<unsigned int> v_out(edges_end.begin(), edges_end.end());
        std::merge(v_in.begin(), v_in.end(),
                   v_out.begin(), v_out.end(),
                   std::inserter(v, v.begin()));
        nof_vertices = v.size();

        // Vertices with edges
        std::set<unsigned int> u_labels_set(u_labels_tmp.begin(), u_labels_tmp.end());
        std::set<unsigned int> v_labels_set(v_labels_tmp.begin(), v_labels_tmp.end());
        u_labels_theo.assign(u_labels_set.begin(), u_labels_set.end());
        v_labels_theo.assign(v_labels_set.begin(), v_labels_set.end());
    }
};
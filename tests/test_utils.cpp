/*!
 * @file
 *
 * @brief Testing some functionalities found in the parameters header.
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

#include "multitensor/parameters.hpp"
#include "multitensor/utils.hpp"
#include "fixtures.hpp"

using namespace multitensor::utils;

BOOST_FIXTURE_TEST_SUITE(tests_utils, fixture_rng)

// Checks the Report class
BOOST_AUTO_TEST_CASE(test_report)
{
    Report results{};
    BOOST_TEST(results.nof_realizations == 0);
    BOOST_TEST(results.vec_L2.size() == 0);
    BOOST_TEST(results.max_L2() == std::numeric_limits<double>::lowest());
    BOOST_TEST(results.vec_iter.size() == 0);
    BOOST_TEST(results.vec_term_reason.size() == 0);
    BOOST_TEST(results.duration == 0);
    BOOST_TEST(results.seed == 0);

    // Likelihood
    double max_L2 = std::numeric_limits<double>::lowest();
    for (size_t i = 0; i < (rng() % 100); i++)
    {
        double rand_L2 = static_cast<double>(rng() % 10000 - 5000);
        results.vec_L2.emplace_back(rand_L2);
        max_L2 = std::max(max_L2, rand_L2);
    }
    BOOST_TEST(max_L2 == results.max_L2());
}

// Checks the default random generator
BOOST_AUTO_TEST_CASE(test_rng)
{
    RandomGenerator double_rng{};
    std::vector<double> vals;
    size_t nof_random_vals = (rng() % 100 + 1);
    for (size_t i = 0; i < nof_random_vals; i++)
    {
        vals.emplace_back(double_rng());
    }
    std::set<double> set_vals(vals.begin(), vals.end());
    BOOST_TEST(set_vals.size() == nof_random_vals);
}

BOOST_AUTO_TEST_SUITE_END()

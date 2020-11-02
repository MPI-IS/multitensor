// Copyright (c) 2019, Max Planck Society / Software Workshop - Max Planck Institute for Intelligent Systems
// Distributed under the GNU GPL license version 3
// See file LICENSE.md or at https://github.com/MPI-IS/multitensor/LICENSE.md

/*!
 * @file
 *
 * @brief Utility classes and functions.
 *
 * @author Jean-Claude Passy (jean-claude.passy@tuebingen.mpg.de)
 * @author Caterina De Bacco (caterina.debacco@tuebingen.mpg.de)
 */

#pragma once

#include <set>
#include <ctime>
#include <random>
#include <boost/random.hpp>

#include "multitensor/tensor.hpp"

namespace multitensor
{

//! Utilities
namespace utils
{

//! @brief Class for containing the results of the algorithm
struct Report
{
    //! The number of realizations
    size_t nof_realizations;

    //! The likelihood of the best results
    std::vector<double> vec_L2;

    //! The number of iterations for each realization
    std::vector<size_t> vec_iter;

    //! The reason for termination
    std::vector<const char *> vec_term_reason;

    //! Duration (in seconds) of the full run
    double duration;

    //! Seed from the random generator
    std::time_t seed;

    //! Default constructor
    Report() = default;

    //! Get maximum likelihood
    double max_L2() const noexcept
    {
        if (vec_L2.size() == 0)
        {
            return std::numeric_limits<double>::lowest();
        }
        else
        {
            auto i = std::max_element(vec_L2.begin(), vec_L2.end());
            return vec_L2[std::distance(vec_L2.begin(), i)];
        }
    }
};

/*!
 * @brief Utility function calculating the number of vertices from input data
 *
 * @tparam vertex_t Type of vertices
 *
 * @param[in] edges_start Labels of vertices where an edge starts
 * @param[in] edges_end Labels of vertices where an edge ends
 *
 * @returns Number of vertices
 */
template <class vertex_t>
size_t get_num_vertices(const std::vector<vertex_t> &edges_start,
                        const std::vector<vertex_t> &edges_end)
{
    std::set<vertex_t> set_e_in(edges_start.begin(), edges_start.end());
    std::set<vertex_t> set_e_out(edges_end.begin(), edges_end.end());
    std::set<vertex_t> set_e;
    std::merge(set_e_in.begin(), set_e_in.end(),
                set_e_out.begin(), set_e_out.end(),
                std::inserter(set_e, set_e.begin()));
    return set_e.size();
}

/*!
 * @brief Class representing a random variate generator
 *
 * @tparam rng_t Type of the random generator
 * @tparam dist_t Type of the distribution
 *
 * @note Built using random number generator together with a random number distribution
 */
template <class rng_t = std::mt19937,
          class dist_t = boost::uniform_real<>>
struct RandomGenerator
{
    //! The seed
    std::time_t seed;

    //! The random number engine;
    rng_t rng;

    //! The random distribution
    dist_t dist;

    /*!
     * @brief Constructor
     *
     * param[in] seed Seed
     */
    RandomGenerator(std::time_t seed = std::time(nullptr))
        : seed(seed)
    {
        rng.seed(static_cast<unsigned int>(seed));
    }

    //! Generate random number
    auto operator()()
    {
        return dist(rng);
    }
};

} // namespace utils
} // namespace multitensor

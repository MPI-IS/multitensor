// Copyright (c) 2019, Max Planck Society / Software Workshop - Max Planck Institute for Intelligent Systems
// Distributed under the GNU GPL license version 3
// See file LICENSE.md or at https://github.com/MPI-IS/multitensor/LICENSE.md

/*!
 * @file
 *
 * @brief Implementation of the solver
 *
 * @author Jean-Claude Passy (jean-claude.passy@tuebingen.mpg.de)
 * @author Caterina De Bacco (caterina.debacco@tuebingen.mpg.de)
 */

#pragma once

#include <iostream>
#include <iomanip>
#include <vector>
#include <cstddef>
#include <memory>

#include "multitensor/graph.hpp"
#include "multitensor/params.hpp"
#include "multitensor/tensor.hpp"
#include "multitensor/utils.hpp"

namespace multitensor
{

//! Solver implementation
namespace solver
{

//! @brief Class for the solver
class Solver
{
private:

    size_t nof_realizations;
    size_t max_nof_iterations;
    size_t nof_convergences;

    /*!
     * @brief Update the u and v components
     *
     * @tparam get_edges_t Proxy type for extracting edges and vertices
     * @tparam affinity_t Affinity matrix type
     * @tparam network_t Network type
     *
     * @param[in] numerator_list List of vertices used to calculate the numerator
     * @param[in] denominator_list List of vertices used to calculate the denominator
     * @param[in] A Network
     * @param[in] w Affinity tensor
     * @param[in] mat_fixed Matrix fixed used for computation
     * @param[in,out] mat_to_update Matrix to update
     * @param[in,out] edges_vertices_proxy Proxy for extracting edges and vertices
     */
    template <class edges_vertices_t,
              class affinity_t,
              class network_t>
    void update_vertices(const std::vector<size_t> &numerator_list,
                         const std::vector<size_t> &denominator_list,
                         const network_t &A,
                         const affinity_t &w,
                         const tensor::Matrix<double> &mat_fixed,
                         tensor::Matrix<double> &mat_to_update,
                         edges_vertices_t edges_vertices_proxy = edges_vertices_t{})
    {
        using direction_t = typename std::decay_t<decltype(A)>::direction_type;

        // Booleans used are compile-time to check if
        //   * we are in the assortative case
        //   * we are using directed network
        constexpr bool assortative = std::is_same_v<affinity_t, tensor::DiagonalTensor<double>> ||
                                     std::is_same_v<affinity_t, tensor::Transpose<tensor::DiagonalTensor<double>>>;
        constexpr bool directed = std::is_same_v<direction_t, boost::bidirectionalS>;

        size_t nof_groups(std::get<0>(w.dims()));
        size_t nof_layers(A.num_layers());
        tensor::Matrix<double> mat_to_update_old(mat_to_update);

        // With the undirected network, references to mat_fixed and mat_to_update are the same
        // so when mat_to_update is updated, mat_fixed will be updated as well.
        // In this case we need to copy of mat_fixed for mat_fixed_old.
        const tensor::Matrix<double> *mat_fixed_old = nullptr;
        if constexpr (directed)
        {
            mat_fixed_old = &mat_fixed;
        }
        else
        {
            mat_fixed_old = new tensor::Matrix<double>(mat_fixed);
        }

        double Z, D, w_k, val_ik, rho_ijkq, Zij_a;

        for (size_t k = 0; k < nof_groups; k++)
        {
            // Calculate Z_u
            Z = 0;
            // Assortative case
            if constexpr (assortative)
            {
                w_k = D = 0;
                for (size_t a = 0; a < nof_layers; a++)
                {
                    w_k += w(k, a);
                }
                for (auto i : denominator_list)
                {
                    D += (*mat_fixed_old)(i, k);
                }
                Z += w_k * D;
            }
            else // Non-assortative case
            {
                for (size_t l = 0; l < nof_groups; l++)
                {
                    w_k = D = 0;
                    for (size_t a = 0; a < nof_layers; a++)
                    {
                        w_k += w(k, l, a);
                    }
                    for (auto i : denominator_list)
                    {
                        D += (*mat_fixed_old)(i, l);
                    }
                    Z += w_k * D;
                }
            }

            if (Z > EPS_PRECISION)
            {
                for (auto i : numerator_list)
                {
                    // Update only if val_ik > 0
                    if (mat_to_update_old(i, k) > EPS_PRECISION)
                    {
                        // Calculate the numerator
                        val_ik = 0;
                        for (size_t a = 0; a < nof_layers; a++)
                        {
                            for (auto its = edges_vertices_proxy.get_edges(i, A(a)); std::get<0>(its) != std::get<1>(its); ++std::get<0>(its))
                            {
                                auto j = edges_vertices_proxy.get_vertex(std::get<0>(its), A(a));

                                // Calculate rho_ijkq
                                rho_ijkq = Zij_a = 0;
                                for (size_t m = 0; m < nof_groups; m++)
                                {
                                    // Assortative case
                                    if constexpr (assortative)
                                    {
                                        Zij_a += mat_to_update_old(i, m) * (*mat_fixed_old)(j, m) * w(m, a);
                                    }
                                    else // Non-assortative case
                                    {
                                        for (size_t l = 0; l < nof_groups; l++)
                                        {
                                            Zij_a += mat_to_update_old(i, m) * (*mat_fixed_old)(j, l) * w(m, l, a);
                                        }
                                    }
                                }
                                if (Zij_a > EPS_PRECISION)
                                {
                                    // Assortative case
                                    if constexpr (assortative)
                                    {
                                        rho_ijkq += (*mat_fixed_old)(j, k) * w(k, a);
                                    }
                                    else // Non-assortative case
                                    {
                                        for (size_t q = 0; q < nof_groups; q++)
                                        {
                                            rho_ijkq += (*mat_fixed_old)(j, q) * w(k, q, a);
                                        }
                                    }
                                    rho_ijkq /= Zij_a;
                                    val_ik += rho_ijkq;
                                }
                            }
                        } // return val_ik

                        mat_to_update(i, k) = mat_to_update_old(i, k) / Z * val_ik;
                        if (std::abs(mat_to_update(i, k)) < EPS_PRECISION)
                        {
                            mat_to_update(i, k) = 0;
                        }
                    }
                }
            }
        } // end cycle over k

        // Free memory if needed
        if constexpr (!directed)
        {
            delete mat_fixed_old;
        }
    }

    /*!
     * @brief Update the affinity matrix
     *
     * @tparam affinity_t Affinity matrix type
     * @tparam network_t Network type
     *
     * @param[in] u_list Indices of vertices with at least one outgoing edge
     * @param[in] v_list Indices of vertices with at least one incoming edge
     * @param[in] A Network
     * @param[in] u Matrix linking vertices in groups for outgoing edges
     * @param[in] v Matrix linking vertices in groups for incoming edges
     * @param[in,out] w Affinity matrix
     *
     * @note Handles the assortative and non-assortative cases.
     * The implementation is a bit redundant but a this single-function implementation
     * has been chosen to maintain consistency with update_vertices
     */
    template <class affinity_t,
              class network_t>
    void update_affinity(const std::vector<size_t> &u_list,
                         const std::vector<size_t> &v_list,
                         const network_t &A,
                         const tensor::Matrix<double> &u,
                         const tensor::Matrix<double> &v,
                         affinity_t &w)
    {
        using graph_t = std::decay_t<decltype(A(0))>;

        size_t nof_groups(std::get<0>(w.dims()));
        size_t nof_vertices(std::get<0>(u.dims()));
        size_t nof_layers(A.num_layers());
        affinity_t w_old(w);

        graph::out_edge_iterator_t<graph_t> eit, eend;
        double Z_kq, Du, Dv, w_kqa, rho_w, Zij_a;

        for (size_t k = 0; k < nof_groups; k++)
        {
            // Assortative case
            if constexpr (std::is_same_v<affinity_t, tensor::DiagonalTensor<double>>)
            {
                // Calculate Z_kq
                Du = Dv = 0;
                for (auto i : v_list)
                {
                    Dv += v(i, k);
                }
                for (auto i : u_list)
                {
                    Du += u(i, k);
                }
                Z_kq = Du * Dv;

                if (Z_kq > EPS_PRECISION)
                {
                    for (size_t a = 0; a < nof_layers; a++)
                    {
                        // Update only if w_ka > 0
                        if (w_old(k, a) > EPS_PRECISION)
                        {

                            w_kqa = 0;
                            for (size_t i = 0; i < nof_vertices; i++)
                            {
                                // Calculate rho_w
                                rho_w = 0;
                                for (tie(eit, eend) = boost::out_edges(i, A(a)); eit != eend; ++eit)
                                {
                                    auto j = boost::target(*eit, A(a)); // OUT-EDGE
                                    Zij_a = 0;
                                    for (size_t m = 0; m < nof_groups; m++)
                                    {
                                        Zij_a += u(i, m) * v(j, m) * w_old(m, a);
                                    }
                                    if (Zij_a > EPS_PRECISION)
                                    {
                                        rho_w += v(j, k) / Zij_a;
                                    }
                                } // return rho_w

                                w_kqa += u(i, k) * rho_w;
                            } // return w_kqa

                            w(k, a) = w_old(k, a) / Z_kq * w_kqa;
                            if (std::abs(w(k, a)) < EPS_PRECISION)
                            {
                                w(k, a) = 0;
                            }
                        }
                    }
                }
            }
            else // Non-assortative case
            {
                for (size_t q = 0; q < nof_groups; q++)
                {
                    // Calculate Z_kq
                    Du = Dv = 0;
                    for (auto i : v_list)
                    {
                        Dv += v(i, q);
                    }
                    for (auto i : u_list)
                    {
                        Du += u(i, k);
                    }
                    Z_kq = Du * Dv;

                    if (Z_kq > EPS_PRECISION)
                    {
                        for (size_t a = 0; a < nof_layers; a++)
                        {
                            // Update only if w_ka > 0
                            if (w_old(k, q, a) > EPS_PRECISION)
                            {

                                w_kqa = 0;
                                for (size_t i = 0; i < nof_vertices; i++)
                                {
                                    // Calculate rho_w
                                    rho_w = 0;
                                    for (tie(eit, eend) = boost::out_edges(i, A(a)); eit != eend; ++eit)
                                    {
                                        auto j = boost::target(*eit, A(a)); // OUT-EDGE
                                        Zij_a = 0;
                                        for (size_t m = 0; m < nof_groups; m++)
                                        {
                                            for (size_t l = 0; l < nof_groups; l++)
                                            {
                                                Zij_a += u(i, m) * v(j, l) * w_old(m, l, a);
                                            }
                                        }
                                        if (Zij_a > EPS_PRECISION)
                                        {
                                            rho_w += v(j, q) / Zij_a;
                                        }
                                    } // return rho_w

                                    w_kqa += u(i, k) * rho_w;
                                } // return w_kqa

                                w(k, q, a) = w_old(k, q, a) / Z_kq * w_kqa;
                                if (std::abs(w(k, q, a)) < EPS_PRECISION)
                                {
                                    w(k, q, a) = 0;
                                }
                            }
                        }
                    }
                } // end cycle over q
            }     // end if assortative
        }         // end cycle over k
    }

    /*!
     * @brief Calculates the likelyhood
     *
     * @tparam affinity_t Affinity matrix type
     * @tparam network_t Network type
     *
     * @param[in] u Matrix linking vertices in groups for outgoing edges
     * @param[in] v Matrix linking vertices in groups for incoming edges
     * @param[in] w Affinity matrix
     * @param[in] A Network
     *
     * @returns Likelyhood
     */
    template <class affinity_t,
              class network_t>
    double calculate_likelyhood(const tensor::Matrix<double> &u,
                                const tensor::Matrix<double> &v,
                                const affinity_t &w,
                                const network_t &A)
    {
        using graph_t = std::decay_t<decltype(A(0))>;

        size_t nof_groups(std::get<0>(w.dims()));
        size_t nof_vertices(std::get<0>(u.dims()));
        size_t nof_layers(A.num_layers());

        double l(0), log_arg, uvw;
        size_t nof_parallel_edges;
        graph::out_edge_iterator_t<graph_t> eit, eend;

        for (size_t alpha = 0; alpha < nof_layers; alpha++)
        {
            for (size_t i = 0; i < nof_vertices; i++)
            {
                for (size_t j = 0; j < nof_vertices; j++)
                {
                    log_arg = 0;
                    for (size_t k = 0; k < nof_groups; k++)
                    {
                        // Assortative case
                        if constexpr (std::is_same_v<affinity_t, tensor::DiagonalTensor<double>>)
                        {
                            uvw = u(i, k) * v(j, k) * w(k, alpha);
                            // Add this term regardeles of the value of A_ijk
                            l -= uvw;
                            // if edge exists, consider this term inside the log argument
                            if (boost::edge(i, j, A(alpha)).second)
                            {
                                log_arg += uvw;
                            }
                        }
                        else
                        {
                            for (size_t q = 0; q < nof_groups; q++)
                            {
                                uvw = u(i, k) * v(j, q) * w(k, q, alpha);
                                // Add this term regardeles of the value of A_ijk
                                l -= uvw;
                                // if edge exists, consider this term inside the log argument
                                if (boost::edge(i, j, A(alpha)).second)
                                {
                                    log_arg += uvw;
                                }
                            }
                        } // end cycle over k and q
                    }

                    if (log_arg > EPS_PRECISION)
                    {
                        nof_parallel_edges = 0;

                        // Cycle over out-neighbors of i in layer a,  --> count parallel edges
                        for (std::tie(eit, eend) = boost::out_edges(i, A(alpha)); eit != eend; ++eit)
                        {
                            if (boost::target(*eit, A(alpha)) == j)
                            {
                                nof_parallel_edges++;
                            }
                        }

                        l += nof_parallel_edges * std::log(log_arg);
                    } // end cycle to count parallel edges
                }     // end cycle over j
            }         // end cycle over i
        }             // end cycle over a
        return l;
    }

    /*!
     * @brief Executes one loop of the solver
     *
     * @tparam affinity_t Affinity matrix type
     * @tparam network_t Network type
     *
     * @param[in] u_list ndices of vertices with at least one outgoing edge
     * @param[in] v_list ndices of vertices with at least one incoming edge
     * @param[in] A Network
     * @param[in,out] u Matrix linking vertices in groups for outgoing edges
     * @param[in,out] v Matrix linking vertices in groups for incoming edges
     * @param[in,out] w Affinity matrix
     * @paran[in,out] iteration Current iteration
     * @paran[in,out] coincide Number of successful checks
     * @paran[in,out] L2 Likelyhoods
     *
     * @returns Whether the algorithm has converged
     */
    template <class affinity_t,
              class network_t>
    termination_reason loop(const std::vector<size_t> &u_list,
                            const std::vector<size_t> &v_list,
                            const network_t &A,
                            tensor::Matrix<double> &u,
                            tensor::Matrix<double> &v,
                            affinity_t &w,
                            size_t &iteration,
                            size_t &coincide,
                            double &L2)
    {
        using direction_t = typename std::decay_t<decltype(A)>::direction_type;

        // Update u (outgoing edges)
        update_vertices<graph::out_edges_target_vertices>(u_list, v_list, A, w, v, u);
        // Update v (incoming edges) if we use directed network
        if constexpr (std::is_same_v<direction_t, boost::bidirectionalS>)
        {

            // Assortative
            if constexpr (std::is_same_v<affinity_t, tensor::DiagonalTensor<double>>)
            {
                update_vertices<graph::in_edges_source_vertices>(v_list, u_list, A, w, u, v);
            }
            else // Non-assortative case
            {
                tensor::Transpose wT(w);
                update_vertices<graph::in_edges_source_vertices>(v_list, u_list, A, wT, u, v);
            }
        }
        // Update w
        update_affinity(u_list, v_list, A, u, v, w);

        // Check for convergence
        if (iteration % 10 == 0)
        {
            double L2_old = L2;
            L2 = calculate_likelyhood(u, v, w, A);
            if (std::abs(L2_old - L2)/std::abs(L2_old) < EPS_PRECISION_LIKELIHOOD)
            {
                coincide++;
            }
            else
            {
                coincide = 0;
            }
        }
        iteration++;

        // Return reason for terminating the loop
        if (coincide == num_conv()) // convergence
        {
            return CONVERGED;
        }
        else if (iteration == max_iter()) // maximum iteration reached
        {
            return MAX_ITER;
        }
        else // otherwise we continue
        {
            return NO_TERMINATION;
        }
    }

public:
    /*!
     * @brief Constructor of a solver
     *
     * @param[in] nof_realizations_ Number of realizations
     * @param[in] max_nof_iterations_ Maximum number of iterations in each realization
     * @param[in] nof_convergences_ Number of successive passed convergence criteria for declaring the results converged
     */
    Solver(size_t nof_realizations, size_t max_nof_iterations, size_t nof_convergences)
        : nof_realizations(nof_realizations),
          max_nof_iterations(max_nof_iterations),
          nof_convergences(nof_convergences)
    {
    }

    /*!
     * @brief Run the solver
     *
     * @tparam affinity_init_t Class type for the initialization of w
     * @tparam affinity_t Class type for w
     * @tparam network_t Network type
     * @tparam random_t Random generator type
     *
     * @param[in] u_list ndices of vertices with at least one outgoing edge
     * @param[in] v_list ndices of vertices with at least one incoming edge
     * @param[in] A Network
     * @param[in] w_init_defined If @c true, w has been initialized
     *            and should not be modified for the first realization
     * @param[in,out] u Matrix linking vertices in groups for outgoing edges
     * @param[in,out] v Matrix linking vertices in groups for incoming edges
     * @param[in,out] w Affinity matrix
     * @param[in,out] random_generator Random generator
     *
     * The solver is an explicit split-operator. At each iteration \f$t\f$, it does the following:
     *   - computes \f$u(t+1) = u^+ = f(u,v,w)\f$
     *   - computes \f$v(t+1) = v^+ = f(u^+,v,w)\f$
     *   - computes \f$w(t+1) = w^+ = f(u^+,v^+,w)\f$
     *
     * Additionally, the check for convergence is performed every 10 iterations
     */
    template <class affinity_init_t,
              class affinity_t,
              class network_t,
              class random_t>
    utils::Report run(const std::vector<size_t> &u_list,
                      const std::vector<size_t> &v_list,
                      const network_t &A,
                      tensor::Matrix<double> &u,
                      tensor::Matrix<double> &v,
                      affinity_t &w,
                      random_t &random_generator,
                      affinity_init_t w_init_obj = affinity_init_t{})
    {
        using direction_t = typename std::decay_t<decltype(A)>::direction_type;

        // Boolean used are compile-time to check if we are using directed network
        constexpr bool directed = std::is_same_v<direction_t, boost::bidirectionalS>;

        // Dimensions
        const size_t nof_groups(std::get<0>(w.dims()));
        const size_t nof_vertices(std::get<0>(u.dims()));
        const size_t nof_layers(A.num_layers());

        // Results
        utils::Report results{};
        results.nof_realizations = num_real();

        // Tensors used within a realization
        tensor::Matrix<double> u_temp(nof_vertices, nof_groups), v_temp;
        affinity_t w_temp(nof_groups, nof_layers), w_init;

        for (size_t i = 0; i < num_real(); i++)
        {
            std::cout << "Running realization # " << i << std::endl;

            // w initialization
            w_init_obj(w, w_temp, random_generator);

            // u and v initialization
            // TO CHECK: there is a difference here in case some nodes in the files dont have in/out edges
            // i.e. remove node 299 in the data file.
            // Compile-time if statement
            // We need to initialize v only if the network is directed
            if constexpr (directed)
            {
                v_temp.resize(nof_vertices, nof_groups);
                initialization::init_tensor_rows_random(v_list, v_temp, random_generator);
            }
            initialization::init_tensor_rows_random(u_list, u_temp, random_generator);

            // Likelihood, convergence criteria and iterations
            double L2(std::numeric_limits<double>::lowest());
            size_t iteration(0), coincide(0);
            assert(max_iter() > 0);

            // Termination reason
            termination_reason term_reason = NO_TERMINATION;
            while (term_reason == NO_TERMINATION)
            {
                // Compile-time if statement
                // Here we choose between directed and undirected network
                if constexpr (directed)
                {
                    term_reason = loop(u_list, v_list, A,
                                       u_temp, v_temp, w_temp,
                                       iteration, coincide, L2);
                }
                else
                {
                    // Undirected network
                    term_reason = loop(u_list, v_list, A,
                                       u_temp, u_temp, w_temp,
                                       iteration, coincide, L2);
                }
            }
            std::cout << "\t... finished after " << iteration << " iterations. "
                      << "Reason: " << get_termination_reason_name(term_reason)
                      << ". Likelihood: " << L2 << std::endl;

            // Update report and best configuration
            results.vec_iter.emplace_back(iteration);
            results.vec_term_reason.emplace_back(get_termination_reason_name(term_reason));
            if (results.max_L2() < L2)
            {
                std::swap(w, w_temp);
                std::swap(u, u_temp);
                if constexpr (directed)
                {
                    std::swap(v, v_temp);
                }
            }
            results.vec_L2.emplace_back(L2);
        }

        // return report
        return results;
    }

    //! @brief Returns the number of realizations.
    size_t num_real() const noexcept
    {
        return nof_realizations;
    }

    //! @brief Returns the number of realizations.
    size_t max_iter() const noexcept
    {
        return max_nof_iterations;
    }

    //! @brief Returns the number convergence checks to declare convergence.
    size_t num_conv() const noexcept
    {
        return nof_convergences;
    }
}; // class Solver

} // namespace solver
} // namespace multitensor

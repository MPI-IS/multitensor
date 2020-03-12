/*!
 * @file
 *
 * @brief Implementation of the solver
 *
 * @author Jean-Claude Passy (jean-claude.passy@tuebingen.mpg.de)
 */

#pragma once

#include <iostream>
#include <iomanip>
#include <vector>
#include <cstddef>

#include "multitensor/graph.hpp"
#include "multitensor/parameters.hpp"
#include "multitensor/tensor.hpp"
#include "multitensor/utils.hpp"

namespace multitensor
{

//! Solver implementation
namespace solver
{

//!@brief Class for the solver
class Solver
{
private:
    unsigned int nof_realizations;
    unsigned int max_nof_iterations;
    unsigned int nof_convergences; // number of times the convergence criterium is satisfied before stopping the simulation

    /*!
     * @brief Update the u component
     *
     * @tparam graph_t Graph type
     *
     * @param[in] u_list Indices of vertices with at least one outgoing edge
     * @param[in] v_list Indices of vertices with at least one incoming edge
     * @param[in] A Network
     * @param[in, out] u Tensor linking vertices in groups for outgoing edges
     * @param[in, out] w_old Affinity tensor (previous step)
     * @param[in, out] u_old Tensor linking vertices in groups for outgoing edges (previous step)
     * @param[in, out] v_old Tensor linking vertices in groups for incoming edges (previous step)
     */
    template <class network_t>
    double update_u(const std::vector<size_t> &u_list,
                    const std::vector<size_t> &v_list,
                    const network_t &A,
                    tensor::Tensor<double> &u,
                    tensor::Tensor<double> &w_old,
                    tensor::Tensor<double> &u_old,
                    tensor::Tensor<double> &v_old)
    {
        using graph_t = std::decay_t<decltype(A(0))>;

        unsigned int nof_groups(std::get<0>(w_old.dims()));
        unsigned int nof_layers(A.num_layers());

        graph::out_edge_iterator<graph_t> eit, eend;
        double Z_u, Du, w_k, dist_u, u_ik, rho_ijkq, Zij_a;
        dist_u = 0;

        for (size_t k = 0; k < nof_groups; k++)
        {
            // Calculate Z_u
            Z_u = 0;
            for (size_t l = 0; l < nof_groups; l++)
            {
                w_k = Du = 0;
                for (size_t a = 0; a < nof_layers; a++)
                {
                    w_k += w_old(k, l, a);
                }
                for (auto i : v_list)
                {
                    Du += v_old(i, l);
                }
                Z_u += w_k * Du;
            }

            if (Z_u > EPS_PRECISION)
            {
                for (auto i : u_list)
                {
                    // Update only if u_ik > 0
                    if (u_old(i, k) > EPS_PRECISION)
                    {
                        // Calculate the numerator
                        u_ik = 0;
                        for (size_t a = 0; a < nof_layers; a++)
                        {
                            for (tie(eit, eend) = boost::out_edges(i, A(a)); eit != eend; ++eit)
                            {
                                graph::Vertex<graph_t> j = target(*eit, A(a)); // OUT-EDGE

                                // Calculate rho_ijkq
                                rho_ijkq = Zij_a = 0;
                                for (size_t m = 0; m < nof_groups; m++)
                                {
                                    for (size_t l = 0; l < nof_groups; l++)
                                    {
                                        Zij_a += u_old(i, m) * v_old(j, l) * w_old(m, l, a);
                                    }
                                }
                                if (Zij_a > EPS_PRECISION)
                                {
                                    for (size_t q = 0; q < nof_groups; q++)
                                    {
                                        rho_ijkq += v_old(j, q) * w_old(k, q, a);
                                    }
                                    rho_ijkq /= Zij_a;
                                }
                                u_ik += rho_ijkq;
                            }
                        } // return u_ik

                        u_ik = u_old(i, k) / Z_u * u_ik;
                        if (std::abs(u_ik) < EPS_PRECISION)
                        {
                            u_ik = 0;
                        }
                        u(i, k) = u_ik;

                        // Calculate max difference
                        dist_u = std::max(std::abs(u(i, k) - u_old(i, k)), dist_u);
                    }
                }
            }
        } // end cycle over k

        // Copy and return
        u_old = u;
        return dist_u;
    }

    /*!
     * @brief Update the v component
     *
     * @tparam graph_t Graph type
     *
     * @param[in] u_list Indices of vertices with at least one outgoing edge
     * @param[in] v_list Indices of vertices with at least one incoming edge
     * @param[in] A Network
     * @param[in, out] v Tensor linking vertices in groups for incoming edges
     * @param[in, out] u Tensor linking vertices in groups for outgoing edges
     * @param[in, out] w_old Affinity tensor (previous step)
     * @param[in, out] v_old Tensor linking vertices in groups for incoming edges (previous step)
     */
    template <class network_t>
    double update_v(const std::vector<size_t> &u_list,
                    const std::vector<size_t> &v_list,
                    const network_t &A,
                    tensor::Tensor<double> &v,
                    tensor::Tensor<double> &u,
                    tensor::Tensor<double> &w_old,
                    tensor::Tensor<double> &v_old)
    {
        using graph_t = std::decay_t<decltype(A(0))>;

        unsigned int nof_groups(std::get<0>(w_old.dims()));
        unsigned int nof_layers(A.num_layers());

        graph::out_edge_iterator<graph_t> eit, eend;
        double Z_v, Dv, w_k, dist_v, v_ik, rho_ijkq, Zij_a;
        dist_v = 0;

        for (size_t k = 0; k < nof_groups; k++)
        {
            // Calculate Z_u
            Z_v = 0;
            for (size_t l = 0; l < nof_groups; l++)
            {
                w_k = Dv = 0;
                for (size_t a = 0; a < nof_layers; a++)
                {
                    w_k += w_old(l, k, a);
                }
                for (auto i : u_list)
                {
                    Dv += u(i, l);
                }
                Z_v += w_k * Dv;
            }

            if (Z_v > EPS_PRECISION)
            {
                for (auto i : v_list)
                {
                    // Update only if u_ik > 0
                    if (v_old(i, k) > EPS_PRECISION)
                    {
                        // Calculate the numerator
                        v_ik = 0;
                        for (size_t a = 0; a < nof_layers; a++)
                        {
                            for (tie(eit, eend) = boost::out_edges(i, A(a)); eit != eend; ++eit)
                            {
                                graph::Vertex<graph_t> j = target(*eit, A(a)); // OUT-EDGE

                                // Calculate rho_ijkq
                                rho_ijkq = Zij_a = 0;
                                for (size_t m = 0; m < nof_groups; m++)
                                {
                                    for (size_t l = 0; l < nof_groups; l++)
                                    {
                                        Zij_a += u(j, m) * v_old(i, l) * w_old(m, l, a);
                                    }
                                }
                                if (Zij_a > EPS_PRECISION)
                                {
                                    for (size_t q = 0; q < nof_groups; q++)
                                    {
                                        rho_ijkq += u(j, q) * w_old(q, k, a);
                                    }
                                    rho_ijkq /= Zij_a;
                                }
                                v_ik += rho_ijkq;
                            }
                        } // return u_ik

                        v_ik = v_old(i, k) / Z_v * v_ik;
                        if (std::abs(v_ik) < EPS_PRECISION)
                        {
                            v_ik = 0;
                        }
                        v(i, k) = v_ik;

                        // Calculate max difference
                        dist_v += std::abs(v(i, k) - v_old(i, k));
                    }
                }
            }
        } // end cycle over k

        // Copy and return
        v_old = v;
        return dist_v;
    }

    /*!
     * @brief Update the affinity matrix
     *
     * @tparam graph_t Graph type
     *
     * @param[in] u_list Indices of vertices with at least one outgoing edge
     * @param[in] v_list Indices of vertices with at least one incoming edge
     * @param[in] A Network
     * @param[in,out] w Affinity tensor
     * @param[in,out] v Tensor linking vertices in groups for incoming edges
     * @param[in,out] u Tensor linking vertices in groups for outgoing edges
     * @param[in,out] w_old Affinity tensor (previous step)
     */
    template <class network_t>
    double update_w(const std::vector<size_t> &u_list,
                    const std::vector<size_t> &v_list,
                    const network_t &A,
                    tensor::Tensor<double> &w,
                    tensor::Tensor<double> &v,
                    tensor::Tensor<double> &u,
                    tensor::Tensor<double> &w_old)
    {
        using graph_t = std::decay_t<decltype(A(0))>;

        unsigned int nof_groups(std::get<0>(w_old.dims()));
        unsigned int nof_vertices(std::get<0>(u.dims()));
        unsigned int nof_layers(A.num_layers());

        graph::out_edge_iterator<graph_t> eit, eend;
        double Z_kq, Du, Dv, dist_w, w_kqa, rho_w, Zij_a;
        dist_w = 0;

        for (size_t k = 0; k < nof_groups; k++)
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
                                    graph::Vertex<graph_t> j = target(*eit, A(a)); // OUT-EDGE
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

                            w_kqa = w_old(k, q, a) / Z_kq * w_kqa;
                            if (std::abs(w_kqa) < EPS_PRECISION)
                            {
                                w_kqa = 0;
                            }
                            w(k, q, a) = w_kqa;

                            // Calculate max difference
                            dist_w = std::max(std::abs(w(k, q, a) - w_old(k, q, a)), dist_w);
                        }
                    }
                }
            } // end cycle over q
        }     // end cycle over k

        // Copy and return
        w_old = w;
        return dist_w;
    }

    /*!
     * @brief Calculates the likelyhood
     *
     * @tparam graph_t Graph type
     *
     * @param[in] u Tensor linking vertices in groups for outgoing edges
     * @param[in] v Tensor linking vertices in groups for incoming edges
     * @param[in] w Affinity tensor
     * @param[in] A Network
     *
     * @returns Likelyhood
     */
    template <class network_t>
    double calculate_likelyhood(const tensor::Tensor<double> &u,
                                const tensor::Tensor<double> &v,
                                const tensor::Tensor<double> &w,
                                const network_t &A)
    {
        using graph_t = std::decay_t<decltype(A(0))>;

        size_t nof_groups(std::get<0>(w.dims()));
        size_t nof_vertices(std::get<0>(u.dims()));
        size_t nof_layers(A.num_layers());

        double l(0), log_arg, uvw;
        unsigned int nof_parallel_edges;
        graph::out_edge_iterator<graph_t> eit, eend;

        for (size_t alpha = 0; alpha < nof_layers; alpha++)
        {
            for (size_t i = 0; i < nof_vertices; i++)
            {
                for (size_t j = 0; j < nof_vertices; j++)
                {
                    log_arg = 0;
                    for (size_t k = 0; k < nof_groups; k++)
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
     * @tparam graph_t Graph type
     *
     * @param[in] u_list ndices of vertices with at least one outgoing edge
     * @param[in] v_list ndices of vertices with at least one incoming edge
     * @param[in] A Network
     * @paran[in] iteration Current iteration
     * @paran[in,out] coincide Number of successful checks
     * @paran[in,out] L2 Likelyhoods
     * @param[in,out] w Affinity tensor
     * @param[in,out] u Tensor linking vertices in groups for outgoing edges
     * @param[in,out] v Tensor linking vertices in groups for incoming edges
     *
     * @returns Whether the algorithm has converged
     */
    template <class network_t>
    termination_reason loop(const std::vector<size_t> &u_list,
                            const std::vector<size_t> &v_list,
                            const network_t &A,
                            unsigned int &iteration,
                            unsigned int &coincide,
                            double &L2,
                            tensor::Tensor<double> &w,
                            tensor::Tensor<double> &u,
                            tensor::Tensor<double> &v)
    {
        using direction_t = typename std::decay_t<decltype(A)>::direction_type;

        // Variables used for copy
        tensor::Tensor<double> w_old(w), u_old(u), v_old(v);

        // Split updates
        update_u(u_list, v_list, A, u, w_old, u_old, v_old);
        if constexpr (std::is_same_v<direction_t, boost::bidirectionalS>)
        {
            update_v(u_list, v_list, A, v, u, w_old, v_old);
        }
        update_w(u_list, v_list, A, w, v, u, w_old);

        // Check for convergence
        if (iteration % 10 == 0)
        {
            double L2_old = L2;
            L2 = calculate_likelyhood(u, v, w, A);
            if (std::abs(L2_old - L2) < EPS_PRECISION)
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
    Solver(unsigned int nof_realizations, unsigned int max_nof_iterations, unsigned int nof_convergences)
        : nof_realizations(nof_realizations),
          max_nof_iterations(max_nof_iterations),
          nof_convergences(nof_convergences)
    {
    }

    /*!
     * @brief Run the solver
     *
     * @tparam w_init_t Class type for the initialization of w
     * @tparam uv_init_t Class type for the initialization of u and v
     * @tparam random_t Random generator type
     * @tparam graph_t Graph type
     *
     * @param[in] u_list ndices of vertices with at least one outgoing edge
     * @param[in] v_list ndices of vertices with at least one incoming edge
     * @param[in] A Network
     * @param[in,out] w Affinity tensor
     * @param[in,out] u Tensor linking vertices in groups for outgoing edges
     * @param[in,out] v Tensor linking vertices in groups for incoming edges
     * @param[in,out] random_generator Random generator
     *
     * The solver is an explicit split-operator. At each iteration \f$t\f$, it does the following:
     *   - computes \f$u(t+1) = u^+ = f(u,v,w)\f$
     *   - computes \f$v(t+1) = v^+ = f(u^+,v,w)\f$
     *   - computes \f$w(t+1) = w^+ = f(u^+,v^+,w)\f$
     *
     * Additionally, the check for convergence is performed every 10 iterations
     */
    template <class w_init_t = initialization::init_tensor_symmetric_random,
              class uv_init_t = initialization::init_tensor_partial_random,
              class random_t,
              class network_t>
    utils::Report run(const std::vector<size_t> &u_list, const std::vector<size_t> &v_list,
                      const network_t &A,
                      tensor::Tensor<double> &w, tensor::Tensor<double> &u, tensor::Tensor<double> &v,
                      random_t &random_generator,
                      w_init_t w_init_obj = w_init_t{},
                      uv_init_t uv_init_obj = uv_init_t{})
    {
        using direction_t = typename std::decay_t<decltype(A)>::direction_type;

        // Dimensions
        const size_t nof_groups(std::get<0>(w.dims()));
        const size_t nof_vertices(std::get<0>(u.dims()));
        const size_t nof_layers(A.num_layers());

        // Results
        utils::Report results{};
        results.nof_realizations = num_real();

        for (unsigned int i = 0; i < num_real(); i++)
        {
            std::cout << "Running realization # " << i << std::endl;

            // Tensors used for this realization
            // u,v,w are used to store the best configuration
            tensor::Tensor<double>
                w_temp(nof_groups, nof_groups, nof_layers),
                u_temp(nof_vertices, nof_groups),
                v_temp;

            // Initialize tensors
            w_init_obj(w_temp, random_generator);
            // TO CHECK: there is a difference here in case some nodes in the files dont have in/out edges
            // i.e. remove node 299 in the data file.
            if constexpr (std::is_same_v<direction_t, boost::bidirectionalS>)
            {
                v_temp.resize(nof_vertices, nof_groups);
                uv_init_obj(v_list, v_temp, random_generator);
            }
            uv_init_obj(u_list, u_temp, random_generator);

            // Likelihood,convergence criteria and iterations
            double L2(std::numeric_limits<double>::max());
            unsigned int iteration(0), coincide(0);
            assert(max_iter() > 0);

            // Termination reason
            termination_reason term_reason = NO_TERMINATION;
            while (term_reason == NO_TERMINATION)
            {
                if constexpr (std::is_same_v<direction_t, boost::bidirectionalS>)
                {
                    term_reason = loop(u_list, v_list, A,
                                       iteration, coincide, L2,
                                       w_temp, u_temp, v_temp);
                }
                else
                {
                    term_reason = loop(u_list, u_list, A,
                                       iteration, coincide, L2,
                                       w_temp, u_temp, v_temp);
                }
            }

            std::cout << "\t... finished after " << iteration << " iterations. "
                      << "Reason: " << get_termination_reason_name(term_reason) << std::endl;

            // Update report and best configuration
            results.vec_iter.emplace_back(iteration);
            results.vec_term_reason.emplace_back(get_termination_reason_name(term_reason));
            if (results.max_L2() < L2)
            {
                std::swap(w, w_temp);
                std::swap(u, u_temp);
                std::swap(v, v_temp);
            }
            results.vec_L2.emplace_back(L2);
        }

        // return report
        return results;
    }

    //! @brief Returns the number of realizations.
    unsigned int num_real() const noexcept
    {
        return nof_realizations;
    }

    //! @brief Returns the number of realizations.
    unsigned int max_iter() const noexcept
    {
        return max_nof_iterations;
    }

    //! @brief Returns the number convergence checks to declare convergence.
    unsigned int num_conv() const noexcept
    {
        return nof_convergences;
    }
};

} // namespace solver
} // namespace multitensor
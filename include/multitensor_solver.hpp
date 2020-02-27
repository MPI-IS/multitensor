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

#include "multitensor_graph.hpp"
#include "multitensor_parameters.hpp"
#include "multitensor_tensor.hpp"

namespace multitensor
{
namespace solver
{

//!@brief Class for the solver
class Solver
{
private:
    unsigned int nof_realizations;
    unsigned int max_nof_iterations;
    unsigned int nof_convergences; // number of times the convergence criterium is satisfied before stopping the simulation

public:
    /*!
     * @brief Constructor of a solver
     *
     * @param[in] nof_realizations_ number of realizations
     */
    Solver(unsigned int nof_realizations_, unsigned int max_nof_iterations_, unsigned int nof_convergences_)
        : nof_realizations(nof_realizations_),
          max_nof_iterations(max_nof_iterations_),
          nof_convergences(nof_convergences_)
    {
    }

    //! @brief Update u
    template <class graph_t>
    double update_u(const std::vector<size_t> &u_list,
                    const std::vector<size_t> &v_list,
                    const std::vector<graph_t> &A,
                    tensor::Tensor<double> &u,
                    tensor::Tensor<double> &w_old,
                    tensor::Tensor<double> &u_old,
                    tensor::Tensor<double> &v_old)
    {
        unsigned int nof_groups(std::get<0>(w_old.dims()));
        unsigned int nof_layers(A.size());

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
            //std::cout << "k = " << k << " Z_u = " << Z_u << std::endl;
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
                            for (tie(eit, eend) = boost::out_edges(i, A[a]); eit != eend; ++eit)
                            {
                                graph::Vertex<graph_t> j = target(*eit, A[a]); // OUT-EDGE

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

    //! @brief Update v
    template <class graph_t>
    double update_v(const std::vector<size_t> &u_list,
                    const std::vector<size_t> &v_list,
                    const std::vector<graph_t> &A,
                    tensor::Tensor<double> &v,
                    tensor::Tensor<double> &u,
                    tensor::Tensor<double> &w_old,
                    tensor::Tensor<double> &u_old,
                    tensor::Tensor<double> &v_old)
    {
        unsigned int nof_groups(std::get<0>(w_old.dims()));
        unsigned int nof_layers(A.size());

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

            //std::cout << "k = " << k << " Z_v = " << Z_v << std::endl;
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
                            for (tie(eit, eend) = boost::out_edges(i, A[a]); eit != eend; ++eit)
                            {
                                graph::Vertex<graph_t> j = target(*eit, A[a]); // OUT-EDGE

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
                        //dist_v = std::max(std::abs(v(i, k) - v_old(i, k)), dist_v);
                        dist_v += std::abs(v(i, k) - v_old(i, k));
                    }
                }
            }
        } // end cycle over k

        // Copy and return
        v_old = v;
        return dist_v;
    }

    //! @brief Update w
    template <class graph_t>
    double update_w(const std::vector<size_t> &u_list,
                    const std::vector<size_t> &v_list,
                    const std::vector<graph_t> &A,
                    tensor::Tensor<double> &w,
                    tensor::Tensor<double> &v,
                    tensor::Tensor<double> &u,
                    tensor::Tensor<double> &w_old)
    {
        unsigned int nof_groups(std::get<0>(w_old.dims()));
        unsigned int nof_nodes(std::get<0>(u.dims()));
        unsigned int nof_layers(A.size());

        graph::out_edge_iterator<graph_t> eit, eend;
        double Z_kq, Du, Dv, w_k, dist_w, w_kqa, rho_w, Zij_a;
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
                            for (size_t i = 0; i < nof_nodes; i++)
                            {
                                // Calculate rho_w
                                rho_w = 0;
                                for (tie(eit, eend) = boost::out_edges(i, A[a]); eit != eend; ++eit)
                                {
                                    graph::Vertex<graph_t> j = target(*eit, A[a]); // OUT-EDGE
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

    //! Calculate likelyhood
    template <class graph_t>
    double calculate_likelyhood(const tensor::Tensor<double> &u,
                                const tensor::Tensor<double> &v,
                                const tensor::Tensor<double> &w,
                                const std::vector<graph_t> &A)
    {

        size_t nof_groups(std::get<0>(w.dims()));
        size_t nof_nodes(std::get<0>(u.dims()));
        size_t nof_layers(A.size());

        double L2(0), l, log_arg, uvw;
        unsigned int nof_parallel_edges;
        graph::out_edge_iterator<graph_t> eit, eend;

        for (size_t alpha = 0; alpha < nof_layers; alpha++)
        {
            for (size_t i = 0; i < nof_nodes; i++)
            {
                for (size_t j = 0; j < nof_nodes; j++)
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
                            if (boost::edge(i, j, A[alpha]).second)
                            {
                                log_arg += uvw;
                            }
                        }
                    } // end cycle over k and q

                    if (log_arg > EPS_PRECISION)
                    {
                        nof_parallel_edges = 0;

                        // Cycle over out-neighbors of i in layer a,  --> count parallel edges
                        for (std::tie(eit, eend) = boost::out_edges(i, A[alpha]); eit != eend; ++eit)
                        {
                            if (boost::target(*eit, A[alpha]) == j)
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

    //! @brief A solver loop
    template <class random_t,
              class graph_t>
    bool loop(const std::vector<size_t> &u_list,
              const std::vector<size_t> &v_list,
              const std::vector<graph_t> &A,
              const unsigned int iteration,
              unsigned int &coincide,
              double &L2,
              tensor::Tensor<double> &w, tensor::Tensor<double> &w_old, tensor::Tensor<double> &w_f,
              tensor::Tensor<double> &u, tensor::Tensor<double> &u_old, tensor::Tensor<double> &u_f,
              tensor::Tensor<double> &v, tensor::Tensor<double> &v_old, tensor::Tensor<double> &v_f,
              random_t &random_generator)
    {

        // Split updates
        double dist_u = update_u(u_list, v_list, A, u, w_old, u_old, v_old);
        double dist_v = update_v(u_list, v_list, A, v, u, w_old, u_old, v_old);
        double dist_w = update_w(u_list, v_list, A, w, v, u, w_old);

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
        return (coincide == nof_convergences);
    }

    //! @brief Run the solver
    template <class random_t,
              class graph_t>
    void run(const size_t &nof_nodes, const size_t &nof_groups, const size_t &nof_layers,
             const std::vector<size_t> &u_list,
             const std::vector<size_t> &v_list,
             const std::vector<graph_t> &A,
             random_t &random_generator)
    {
        // Affinity tensors for the groups
        tensor::Tensor<double> w(nof_groups, nof_groups, nof_layers);
        tensor::Tensor<double> w_old(nof_groups, nof_groups, nof_layers);
        tensor::Tensor<double> w_f(nof_groups, nof_groups, nof_layers);

        // Matrices linking vertices in groups
        tensor::Tensor<double> u(nof_nodes, nof_groups), v(nof_nodes, nof_groups);
        tensor::Tensor<double> u_old(nof_nodes, nof_groups), v_old(nof_nodes, nof_groups);
        tensor::Tensor<double> u_f(nof_nodes, nof_groups), v_f(nof_nodes, nof_groups);

        for (unsigned int i = 0; i < nof_realizations; i++)
        {
            std::cout << "Running realization " << i << std::endl;

            // Initialize tensors
            w.randomize_symmetric(random_generator);

            // TO CHECK: there is a difference here in case some nodes in the files dont have in/out edges
            // i.e. remove node 299 in the data file.
            v.randomize_partial(random_generator, v_list);
            u.randomize_partial(random_generator, u_list);

            // Update old variables
            u_old = u;
            v_old = v;
            w_old = w;

            // Convergence criteria and iterations
            bool convergence = false;
            unsigned int iteration(0), coincide(0);
            double L2 = std::numeric_limits<double>::lowest();

            while (!convergence & (iteration < max_nof_iterations))
            {
                convergence = loop(u_list, v_list, A,
                                   iteration, coincide, L2,
                                   w, w_old, w_f,
                                   u, u_old, u_f,
                                   v, v_old, v_f,
                                   random_generator);
                /*
                std::cout << std::setprecision(12)
                          << "Iteration " << iteration
                          << ", L2 = " << L2 << std::endl;
                */
                iteration++;
            }

            std::cout << std::setprecision(12)
                      << "Final Likelihood = " << L2
                      << ", number of  iterations:" << iteration << std::endl;

            // Output results
            // w
            std::cout << "Final affinity matrix:" << std::endl;
            for (size_t alpha = 0; alpha < nof_layers; alpha++)
            {
                for (size_t k = 0; k < nof_groups; k++)
                {
                    for (size_t q = 0; q < nof_groups; q++)
                    {
                        std::cout << w(k, q, alpha) << " ";
                    }
                    std::cout << std::endl;
                }
                std::cout << std::endl;
            }
        }
    }

    //! @brief Returns the number of realizations.
    unsigned int num_real() const noexcept
    {
        return nof_realizations;
    }

    //! @brief Returns the number of realizations.
    unsigned int max_iterations() const noexcept
    {
        return max_nof_iterations;
    }
}; // namespace solver
} // namespace solver
} // namespace multitensor
// Copyright (c) 2019, Max Planck Society / Software Workshop - Max Planck Institute for Intelligent Systems
// Distributed under the GNU GPL license version 3
// See file LICENSE.md or at https://github.com/MPI-IS/multitensor/LICENSE.md

/*!
 * @file
 *
 * @brief Some global parameters for the projet.
 *
 * @author Jean-Claude Passy (jean-claude.passy@tuebingen.mpg.de)
 * @author Caterina De Bacco (caterina.debacco@tuebingen.mpg.de)
 */

#pragma once

#include <cstddef>

//! @brief Stringification macros
#define stringify(x) #x
#define xstringify(x) stringify(x)

#ifdef PROJECT_VERSION
//! @brief Version number
#define MULTITENSOR_VERSION xstringify(PROJECT_VERSION)
#else
#define MULTITENSOR_VERSION "Unknown"
#endif

//! @brief Type for the dimension
typedef size_t dimension_t;

//! @brief Absolute precision
const double EPS_PRECISION = 1e-6;

//! @brief Absolute precision for logLikelihood convergence
const double EPS_PRECISION_LIKELIHOOD = 1e-4;

//! @brief Noise magnitude
const double EPS_NOISE = 0.1;

//! @brief Enumeration of the reasons for termination
enum termination_reason
{
    NO_TERMINATION,
    MAX_ITER,
    CONVERGED,
};

//! @brief Helper to display termination reason
inline const char *
get_termination_reason_name(termination_reason tr)
{
    switch (tr)
    {
    case MAX_ITER:
        return "MAX_ITER";
    case CONVERGED:
        return "CONVERGED";
    case NO_TERMINATION:
    default:
        return "NO_TERMINATION";
    }
}

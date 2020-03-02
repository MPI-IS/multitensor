/*!
 * @file
 *
 * @brief Some global parameters for the projet.
 *
 * @author Jean-Claude Passy (jean-claude.passy@tuebingen.mpg.de)
 */

#pragma once

// Stringification macros
#define stringify(x) #x
#define xstringify(x) stringify(x)

//! @brief Version number
#ifdef PROJECT_VERSION
#define MULTITENSOR_VERSION xstringify(PROJECT_VERSION)
#else
#define MULTITENSOR_VERSION "Unknown"
#endif

//! @brief Type for the dimension
typedef size_t dimension_t;

//! @brief Absolute precision
const double EPS_PRECISION = 1e-6;

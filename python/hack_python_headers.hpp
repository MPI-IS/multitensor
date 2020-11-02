// Copyright (c) 2019, Max Planck Society / Software Workshop - Max Planck Institute for Intelligent Systems
// Distributed under the GNU GPL license version 3
// See file LICENSE.md or at https://github.com/MPI-IS/multitensor/LICENSE.md

/*!
 * @file
 *
 * @brief File for hacking the Python headers
 *
 * @author Jean-Claude Passy (jean-claude.passy@tuebingen.mpg.de)
 */

#ifndef MULTITENSOR_HACK_PYTHON_H__
#define MULTITENSOR_HACK_PYTHON_H__

#if defined(_WIN32) && defined(_DEBUG)
#define _DEBUG_IS_DEFINED
#undef _DEBUG
#endif

#include <Python.h>

#if defined(_DEBUG_IS_DEFINED)
#define _DEBUG
#undef _DEBUG_IS_DEFINED
#endif

#endif /* MULTITENSOR_HACK_PYTHON_H__ */

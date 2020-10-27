/*!
 * @file
 *
 * @brief File for hacking the Python headers
 *
 * @author Jean-Claude Passy (jean-claude.passy@tuebingen.mpg.de)
 *
 * This file modifying the inclusion of the python libs on Windows
 * and prevents the linking to the debug libraries python3x_d.lib
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

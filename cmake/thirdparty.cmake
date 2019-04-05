# Copyright 2019, Max Planck Society.
# Distributed under the "GNU GPL v3" licence.
# (See accompanying file LICENSE.md)

# Boost
# disable auto link
add_definitions(-DBOOST_ALL_NO_LIB)
if(NOT DEFINED Boost_USE_STATIC_LIBS)
    set(Boost_USE_STATIC_LIBS ON) # More convenient
endif()

if(NOT Boost_USE_STATIC_LIBS)
    # link against dynamic libraries
    add_definitions(-DBOOST_ALL_DYN_LINK)
endif()


set(Boost_REALPATH              ON)
set(Boost_USE_MULTITHREADED     ON)
set(Boost_DEBUG                 ON)
set(Boost_DETAILED_FAILURE_MSG  ON)
set(Boost_NO_BOOST_CMAKE        ON)

set(Boost_ADDITIONAL_VERSIONS
    "1.54" "1.54.0" "1.55" "1.55.0" "1.56" "1.56.0" "1.57.0" "1.58" "1.58.0" "1.59" "1.59.0"
    "1.60" "1.60.0" "1.61" "1.61.0" "1.62" "1.62.0" "1.63" "1.63.0" "1.64" "1.64.0" "1.65" "1.65.0" "1.68" "1.68.0" "1.69" "1.69.0")

# if BOOST_ROOT is defined, we look only there, otherwise we also look into the system paths
if(DEFINED BOOST_ROOT)
    set(Boost_NO_SYSTEM_PATHS ON)
else()
    set(Boost_NO_SYSTEM_PATHS OFF)
endif()

find_package(Boost COMPONENTS unit_test_framework system thread program_options system filesystem)

if(NOT ${Boost_FOUND})
    message(FATAL_ERROR "[BOOST] Boost not found. Please set BOOST_ROOT in your command line.")
endif()


# Doxygen
find_package(Doxygen)


# Python3
set(MULTI_TENSOR_PYTHON_EXTENSIONS FALSE)
if(ENABLE_PYTHON_WRAPPER)
    message(STATUS "[PYTHON] Configuring Python bindings")

    # Find python3
    set (CMAKE_FIND_FRAMEWORK NEVER) # needed on OSX to point the the virtualenv
    find_package(Python3 COMPONENTS Interpreter)
    if(Python3_FOUND)
        message(STATUS "[PYTHON] python3 interpreter found at following location '${Python3_EXECUTABLE}'")
        message(STATUS "[PYTHON] python3 libraries found, include '${Python3_INCLUDE_DIRS}'")

        # Determining the Python shared library extensions
        if(NOT DEFINED PYTHON_SHARED_LIBRARY_EXTENSION)
            execute_process(
                COMMAND ${Python3_EXECUTABLE} -c "import imp; print([i[0] for i in imp.get_suffixes() if i[2] == imp.C_EXTENSION][0])"
                OUTPUT_VARIABLE PYTHON_SHARED_LIBRARY_EXTENSION
                ERROR_VARIABLE  PYTHON_SHARED_LIBRARY_EXTENSION_ERR
                OUTPUT_STRIP_TRAILING_WHITESPACE
            )

            if((PYTHON_SHARED_LIBRARY_EXTENSION_ERR STREQUAL "") AND NOT ("${PYTHON_SHARED_LIBRARY_EXTENSION}" STREQUAL ""))
                set(PYTHON_SHARED_LIBRARY_EXTENSION ${PYTHON_SHARED_LIBRARY_EXTENSION} CACHE STRING "Python shared library extensions" FORCE)
            else()
                unset(PYTHON_SHARED_LIBRARY_EXTENSION)
                message(WARNING "Cannot determine python3 extension suffix")
            endif()
        endif()

        #[[ For now we do not need neither numpy nor cython
        if(NOT DEFINED NUMPY_INCLUDE_PATH)
            # Determining the numpy headers location
            execute_process(
                COMMAND ${Python3_EXECUTABLE} -c "import numpy; print(numpy.get_include())"
                OUTPUT_VARIABLE NUMPY_INCLUDE_PATH
                ERROR_VARIABLE  NUMPY_ERROR
                OUTPUT_STRIP_TRAILING_WHITESPACE
            )

            if((NUMPY_ERROR STREQUAL "") AND NOT ("${NUMPY_INCLUDE_PATH}" STREQUAL "") AND EXISTS "${NUMPY_INCLUDE_PATH}")
                set(NUMPY_INCLUDE_PATH "${NUMPY_INCLUDE_PATH}" CACHE STRING "Python numpy include folder" FORCE)
            else()
                unset(NUMPY_INCLUDE_PATH)
                message(WARNING "Numpy header location not found: is numpy installed? include: '${NUMPY_INCLUDE_PATH}' / error: '${NUMPY_ERROR}'")
            endif()
        endif()

        # Find Cython
        find_program(CYTHON_PROGRAM NAMES cython)
        if(CYTHON_PROGRAM)
            message(STATUS "[PYTHON] Cython wrapper found ${CYTHON_PROGRAM}")
        else()
            message(STATUS "[PYTHON] Cython wrapper not found - will use the .cpp file provided")
        endif()
        ]]

    endif()

endif()

# if all the conditions are met, we enable the python+cython+numpy extensions
if((NOT "${PYTHON_SHARED_LIBRARY_EXTENSION}" STREQUAL ""))
    set(MULTI_TENSOR_PYTHON_EXTENSIONS TRUE)
endif()

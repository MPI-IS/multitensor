# Copyright (c) 2019, Max Planck Society / Software Workshop - Max Planck Institute for Intelligent Systems
# Distributed under the GNU GPL license version 3
# See file LICENSE.md or at https://github.com/MPI-IS/multitensor/LICENSE.md
#
# This file contains the logic for setting up the thridparties configurations
#

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
    "1.60" "1.60.0" "1.61" "1.61.0" "1.62" "1.62.0" "1.63" "1.63.0" "1.64" "1.64.0" "1.65" "1.65.0" "1.68" "1.68.0"
    "1.69" "1.69.0" "1.70" "1.70.0" "1.71" "1.71.0" "1.72" "1.72.0" "1.73" "1.73.0" "1.74" "1.75.0")

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
if(NOT DOXYGEN_FOUND)
    message(STATUS "[DOC] Doxygen not found")
endif()

# Python3 extension
set(MULTITENSOR_PYTHON_EXTENSIONS FALSE)
set(MULTITENSOR_PYTHON_SPHINX FALSE)
set(MULTITENSOR_PYTHON_TEST_NOSE FALSE)

if(ENABLE_PYTHON_WRAPPER)
    message(STATUS "[PYTHON] Checking python3 configuration")

    # Find python3
    set(CMAKE_FIND_FRAMEWORK NEVER) # needed on OSX to point to the virtualenv
    find_package (Python3 COMPONENTS Interpreter Development)
    if(Python3_FOUND)
        message(STATUS "[PYTHON] python3 interpreter found at following location '${Python3_EXECUTABLE}'")
        message(STATUS "[PYTHON] python3 libraries found, include '${Python3_INCLUDE_DIRS}'")

        # Determining the Python shared library extensions
        if(NOT DEFINED PYTHON_SHARED_LIBRARY_EXTENSION)
            execute_process(
                COMMAND ${Python3_EXECUTABLE} -c "import importlib.machinery as mach; print(mach.EXTENSION_SUFFIXES[0])"
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
        get_filename_component(_python_path "${Python3_EXECUTABLE}" DIRECTORY)
        find_program(CYTHON_PROGRAM
            NAMES cython
            HINTS ${_python_path})
        if(CYTHON_PROGRAM)
            message(STATUS "[PYTHON] Cython wrapper found ${CYTHON_PROGRAM}")
        else()
            message(STATUS "[PYTHON] Cython wrapper not found - will use the .cpp file provided")
        endif()
    else()
      message(STATUS "[PYTHON] Python not found!")
    endif()

    # if all the conditions are met, we enable the python+cython+numpy extensions
    if((NOT "${NUMPY_INCLUDE_PATH}" STREQUAL "") AND (NOT "${PYTHON_SHARED_LIBRARY_EXTENSION}" STREQUAL ""))
        set(MULTITENSOR_PYTHON_EXTENSIONS TRUE)

        get_filename_component(_python_path "${Python3_EXECUTABLE}" DIRECTORY)

        # Sphinx documentation
        find_program(SPHINX_EXECUTABLE
            NAMES sphinx-build
            PATHS ${_python_path})
        if(SPHINX_EXECUTABLE)
            message(STATUS "[PYTHON] Sphinx found ${SPHINX_EXECUTABLE}")
            set(MULTITENSOR_PYTHON_SPHINX TRUE)
        else()
            message(STATUS "[PYTHON] Sphinx not found - cannot build Sphinx documentation")
        endif()


        # If possible, use nosetests to get an XML report.
        # Find nose
        find_program(NOSE_PROGRAM
            NAMES nosetests
            HINTS ${_python_path})
        if(NOSE_PROGRAM)
            message(STATUS "[PYTHON] Nosetests found ${NOSE_PROGRAM}")
            set(MULTITENSOR_PYTHON_TEST_NOSE TRUE)
        else()
            message(STATUS "[PYTHON] Nosetests not found - using standard unittest")
        endif()
    endif()

endif()


# Benchmark
set(MULTITENSOR_ENABLE_BENCHMARK FALSE)
if(ENABLE_BENCHMARK)
    message(STATUS "[BENCHMARK] Checking python3 configuration")

    # Python might have been already found for the extension
    if(NOT Python3_Interpreter_FOUND)
        set(CMAKE_FIND_FRAMEWORK NEVER) # needed on OSX to point to the virtualenv
        find_package (Python3 COMPONENTS Interpreter)
    endif()
    if(Python3_Interpreter_FOUND)
        message(STATUS "[BENCHMARK] python3 interpreter found at following location '${Python3_EXECUTABLE}'")
        set(MULTITENSOR_ENABLE_BENCHMARK TRUE)
    else()
        message(STATUS "[BENCHMARK] python3 interpreter not found")
    endif()
endif()

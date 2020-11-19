# Copyright (c) 2019, Max Planck Society / Software Workshop - Max Planck Institute for Intelligent Systems
# Distributed under the GNU GPL license version 3
# See file LICENSE.md or at https://github.com/MPI-IS/multitensor/LICENSE.md
#
# This file contains the logic for building the python extension and running the tests
#

if(NOT MULTITENSOR_PYTHON_EXTENSIONS)
    return()
endif()

message(STATUS "[PYTHON] Configuring Python bindings")

# Paths
set(PYTHON_SRC_DIR ${CMAKE_SOURCE_DIR}/python)
set(PYTHON_PACKAGE_SRC_DIR ${PYTHON_SRC_DIR}/package)

# either generate the .cpp from the .pyx file is Cython exists
# or copy the .cpp from the repo
set(_cython_generated_file_bare "${CMAKE_CURRENT_BINARY_DIR}/multitensor_bare.cpp")
set(_cython_generated_file "${CMAKE_CURRENT_BINARY_DIR}/multitensor.cpp")
if(CYTHON_PROGRAM)
    add_custom_command(
        OUTPUT
            ${_cython_generated_file_bare}
        COMMAND
            ${CYTHON_PROGRAM}
        ARGS
            --cplus -3
            "${PYTHON_PACKAGE_SRC_DIR}/multitensor.pyx"
            -o
            ${_cython_generated_file_bare}
        COMMENT
            "[PYTHON] Generating the binding files with Cython"
        DEPENDS
            "${PYTHON_PACKAGE_SRC_DIR}/multitensor.pyx"
            ${multitensor_src}
    )
else()
    add_custom_command(
        OUTPUT
            ${_cython_generated_file_bare}
        COMMAND
            ${CMAKE_COMMAND}
        ARGS
            -E copy
            ${PYTHON_PACKAGE_SRC_DIR}/multitensor.cpp
            ${_cython_generated_file_bare}
        COMMENT
            "[PYTHON] Copying the .ccp file form the repository"
        DEPENDS
            ${PYTHON_PACKAGE_SRC_DIR}/multitensor.cpp
    )
endif()

# Hack the python headers for builds on Windows
if(MSVC)
    add_custom_command(
        OUTPUT
            ${_cython_generated_file}
        COMMAND ${CMAKE_COMMAND}
            ARGS
                "-Dfile1=${PYTHON_SRC_DIR}/hack_python_headers.hpp"
                "-Dfile2=${_cython_generated_file_bare}"
                "-Doutput=${_cython_generated_file}"
                -P "${CMAKE_SOURCE_DIR}/cmake/cat_files.cmake"
        COMMENT "[PYTHON] Windows builds - hacking Python headers"
        DEPENDS
            ${PYTHON_SRC_DIR}/hack_python_headers.hpp
            ${_cython_generated_file_bare})
else()
    add_custom_command(
        OUTPUT
            ${_cython_generated_file}
        COMMAND ${CMAKE_COMMAND}
            ARGS
                -E copy
                ${_cython_generated_file_bare}
                ${_cython_generated_file}
        COMMENT "[PYTHON] Non-Windows builds - copying cpp file"
        DEPENDS
            ${_cython_generated_file_bare})
endif()

add_custom_target(
    cython_pyx_file
    SOURCES "${PYTHON_PACKAGE_SRC_DIR}/multitensor.pyx")
set_target_properties(
    cython_pyx_file
    PROPERTIES
        FOLDER "Python")

# Shared library
add_library(multitensor_py
    SHARED ${_cython_generated_file})
target_include_directories(
    multitensor_py
    PRIVATE
        ${NUMPY_INCLUDE_PATH}
        ${Python3_INCLUDE_DIR})
target_compile_definitions(multitensor_py
    PRIVATE -DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)
target_link_libraries(multitensor_py
    PRIVATE multitensor multitensor_utils Python3::Python Boost::filesystem)
set_target_properties(multitensor_py
    PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        OUTPUT_NAME "multitensor"
        PREFIX ""
        SUFFIX "${PYTHON_SHARED_LIBRARY_EXTENSION}"
        FOLDER "Python")

# If macOS, add MACOSX_DEPLOYMENT_TARGET
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")

    execute_process(
        COMMAND ${Python3_EXECUTABLE} -c
            "import sysconfig; print(sysconfig.get_config_var('MACOSX_DEPLOYMENT_TARGET'))"
        OUTPUT_VARIABLE MACOSX_TARGET
        ERROR_VARIABLE  MACOSX_TARGET_ERROR
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(MACOSX_TARGET_ERROR STREQUAL "")
        target_compile_options(multitensor_py PUBLIC "-mmacosx-version-min=${MACOSX_TARGET}")
        target_link_options(multitensor_py PUBLIC "-mmacosx-version-min=${MACOSX_TARGET}")
    else()
        message(WARNING "[PYTHON] MACOSX_DEPLOYMENT_TARGET not found")
    endif()
endif()

# Setup.py to build distributions
set(_setup_py_in ${PYTHON_SRC_DIR}/setup.py)

add_custom_target(multitensor_py_dist ALL
    DEPENDS ${_setup_py_in})
add_custom_command(
    TARGET multitensor_py_dist
    POST_BUILD
    WORKING_DIRECTORY $<TARGET_FILE_DIR:multitensor_py>
    COMMAND ${CMAKE_COMMAND} -E copy ${_setup_py_in} $<TARGET_FILE_DIR:multitensor_py>
    COMMAND ${Python3_EXECUTABLE} ${PYTHON_SRC_DIR}/update_setup.py setup.py
    ${CMAKE_PROJECT_VERSION} $<TARGET_FILE_NAME:multitensor_py>)
add_dependencies(multitensor_py_dist multitensor_py)

# Functional tests.
# Copy file
file(COPY ${PYTHON_SRC_DIR}/tests.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
# If possible, run with nose to get XML report
if(MULTITENSOR_PYTHON_TEST_NOSE)
    add_test(
        NAME multitensor_python_test
        COMMAND ${NOSE_PROGRAM} --with-xunit -v
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
else()
    add_test(
        NAME multitensor_python_test
        COMMAND ${Python3_EXECUTABLE} -m unittest -vvv
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif()

# Add path to the library
set_tests_properties(multitensor_python_test
    PROPERTIES
    ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:multitensor_py>"
    ENVIRONMENT "ROOT_DIR=${CMAKE_SOURCE_DIR}")

# Sphinx documentation for Python extension
if(MULTITENSOR_PYTHON_SPHINX)
    set(SPHINX_SRC_DIR ${PYTHON_SRC_DIR}/sphinx_doc)
    set(SPHINX_TARGET_DIR ${CMAKE_CURRENT_BINARY_DIR}/sphinx_doc)
    file(COPY ${SPHINX_SRC_DIR} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

    # To get version and project name from cmake
    configure_file(
        ${SPHINX_SRC_DIR}/source/conf.py
        ${SPHINX_TARGET_DIR}/source/conf.py)

    add_custom_target(sphinx
        COMMAND
            ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/sphinx_doc/build
        COMMAND
        PYTHONPATH=$<TARGET_FILE_DIR:multitensor_py>
            ${SPHINX_EXECUTABLE} -q -b html
            ${SPHINX_TARGET_DIR}/source
            ${SPHINX_TARGET_DIR}/build
        COMMENT "[PYTHON] Generating Sphinx documentation"
        DEPENDS multitensor_py
        SOURCES
            ${SPHINX_SRC_DIR}/source/conf.py
            ${SPHINX_SRC_DIR}/source/definitions.rst
            ${SPHINX_SRC_DIR}/source/index.rst
    )
    set_target_properties(sphinx
        PROPERTIES
            FOLDER "Documentation")
endif()

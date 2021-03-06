# Copyright (c) 2019, Max Planck Society / Software Workshop - Max Planck Institute for Intelligent Systems
# Distributed under the GNU GPL license version 3
# See file LICENSE.md or at https://github.com/MPI-IS/multitensor/LICENSE.md


cmake_minimum_required (VERSION 3.5)
cmake_policy(SET CMP0054 NEW)

if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.13")
    cmake_policy(SET CMP0076 OLD)
endif()

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Extract version number
file(READ ${CMAKE_SOURCE_DIR}/version.txt ver)
project(MultiTensor VERSION ${ver})
message(STATUS "Project ${CMAKE_PROJECT_NAME} v${CMAKE_PROJECT_VERSION}")

enable_testing()

option(ENABLE_PYTHON_WRAPPER "Enable the python extensions" ON)
option(ENABLE_BENCHMARK "Enable the benchmark" ON)

# 0:silent
add_definitions(-DMULTITENSOR_VERBOSE=1)

#################################
# thirdparties and configurations
include(cmake/compiler_settings.cmake)
include(cmake/thirdparty.cmake)


################################
# Main libraries and executables
set(multitensor_src
    ${CMAKE_SOURCE_DIR}/include/multitensor/graph.hpp
    ${CMAKE_SOURCE_DIR}/include/multitensor/initialization.hpp
    ${CMAKE_SOURCE_DIR}/include/multitensor/params.hpp
    ${CMAKE_SOURCE_DIR}/include/multitensor/solver.hpp
    ${CMAKE_SOURCE_DIR}/include/multitensor/tensor.hpp
    ${CMAKE_SOURCE_DIR}/include/multitensor/utils.hpp
    ${CMAKE_SOURCE_DIR}/include/multitensor/main.hpp)

# libraries
add_library(multitensor INTERFACE)
target_include_directories(multitensor
    INTERFACE
        ${CMAKE_SOURCE_DIR}/include/
        ${CMAKE_SOURCE_DIR}/src/)
target_sources(multitensor
    INTERFACE
        $<BUILD_INTERFACE:${multitensor_src} >) # for the IDE to see the library

add_library(multitensor_utils
    ${CMAKE_SOURCE_DIR}/applications/include/app_params.hpp
    ${CMAKE_SOURCE_DIR}/applications/include/app_utils.hpp
    ${CMAKE_SOURCE_DIR}/applications/src/app_utils.cpp)
target_link_libraries(multitensor_utils PUBLIC multitensor Boost::filesystem)
target_include_directories(multitensor_utils
    PUBLIC
    ${CMAKE_SOURCE_DIR}/applications/include/
    ${CMAKE_SOURCE_DIR}/applications/src/)

# executables
add_executable(Multitensor applications/src/multitensor.cpp)
target_link_libraries(Multitensor multitensor multitensor_utils)

# Unit tests
add_executable(test_multitensor
    tests/fixtures.hpp
    tests/test_graph.cpp
    tests/test_initialization.cpp
    tests/test_multitensor.cpp
    tests/test_solver.cpp
    tests/test_tensor.cpp
    tests/test_utils.cpp)
target_link_libraries(test_multitensor
    multitensor
    Boost::unit_test_framework
    Boost::system
    Boost::filesystem)
add_test(
    NAME multitensor_unit_tests
    COMMAND test_multitensor)

# Functional tests
# Helper for command line run
function(add_test_cmd_line
    test_name
    test_program)

    # Set test folder and copy adjacency file
    set(test_folder ${CMAKE_BINARY_DIR}/tests/${test_name})
    file(COPY ${CMAKE_SOURCE_DIR}/data/${test_name}/ DESTINATION ${test_folder}/)
    file(COPY ${CMAKE_SOURCE_DIR}/tests/data/compare.sh DESTINATION ${test_folder}/)
    file(COPY ${CMAKE_SOURCE_DIR}/tests/data/${test_name}/ DESTINATION ${test_folder}/)

    # Run
    add_test(
        NAME "command-line-${test_name}-run"
        COMMAND ${test_program} ${ARGN}
        WORKING_DIRECTORY ${test_folder}
    )

    if("${test_program}" STREQUAL "Multitensor")
        add_test(
            NAME "command-line-${test_name}-verify-content"
            COMMAND ${CMAKE_COMMAND}
                "-Dresults=${test_folder}/results"
                "-Dground_truth=${CMAKE_SOURCE_DIR}/tests/data/${test_name}"
                -P ${CMAKE_SOURCE_DIR}/cmake/verify_content.cmake
        )

        set_tests_properties(
            "command-line-${test_name}-verify-content"
            PROPERTIES
            DEPENDS "command-line-${test_name}-run"
        )
    endif()

endfunction()

add_test_cmd_line(
    "main"
    Multitensor
    --a adjacency.dat
    --k 2
    --s 5489
    --maxit 100)

add_test_cmd_line(
    "undirected"
    Multitensor
    --a adjacency.dat
    --k 2
    --s 5489
    --undirected
    --maxit 100)

add_test_cmd_line(
    "w_input"
    Multitensor
    --a adjacency_k2L4.dat
    --k 2
    --s 5489
    --w w_k2_k2L4_A.dat
    --maxit 100)

add_test_cmd_line(
    "multi_real"
    Multitensor
    --a adjacency_k2L4.dat
    --k 2
    --s 5489
    --w w_k2_k2L4_r2.dat
    --r 2
    --maxit 100)

add_test_cmd_line(
    "assortative"
    Multitensor
    --a adjacency_assortative_k3L4.dat
    --k 3
    --s 5489
    --assortative
    --maxit 100)


#########################
# Python Bindings + Tests
include(cmake/python.cmake)


#########################
# Documentation (Doxygen)
include(cmake/doxygen_doc.cmake)


###########
# Benchmark
include(cmake/benchmark.cmake)

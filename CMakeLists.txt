# Copyright 2019, Max Planck Society.
# Distributed under the "GNU GPL v3" licence.
# (See accompanying file LICENSE.md)

cmake_minimum_required (VERSION 3.5)
cmake_policy(SET CMP0054 NEW)

if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.13")
    cmake_policy(SET CMP0076 OLD)
endif()

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

project(MultiTensor)

enable_testing()

option(ENABLE_PYTHON_WRAPPER "Enable the python extensions" ON)

# 0:silent
add_definitions(-DMULTITENSOR_VERBOSE=1)

#################################
# thirdparties and configurations
include(cmake/compiler_settings.cmake)
include(cmake/thirdparty.cmake)


################################
# Main libraries and executables
set(multitensor_src
    ${CMAKE_SOURCE_DIR}/include/multitensor_graph.hpp
    ${CMAKE_SOURCE_DIR}/include/multitensor_parameters.hpp
    ${CMAKE_SOURCE_DIR}/include/multitensor_solver.hpp
    ${CMAKE_SOURCE_DIR}/include/multitensor_tensor.hpp
    ${CMAKE_SOURCE_DIR}/include/multitensor.hpp)
set(multitensor_legacy_src
    ${CMAKE_SOURCE_DIR}/legacy/include/main.hpp
    ${CMAKE_SOURCE_DIR}/legacy/include/mlg.hpp
    ${CMAKE_SOURCE_DIR}/legacy/src/tools_assortative.cpp
    ${CMAKE_SOURCE_DIR}/legacy/src/tools.cpp
    )
set(multitensor_legacy_default_src
    ${CMAKE_SOURCE_DIR}/legacy/include/cycle_over_realizations.hpp
    ${CMAKE_SOURCE_DIR}/legacy/src/cycle_over_realizations.cpp
    ${CMAKE_SOURCE_DIR}/legacy/src/em_update.cpp
    )
set(multitensor_legacy_undirected_src
    ${CMAKE_SOURCE_DIR}/legacy/include/cycle_over_realizations_undirected.hpp
    ${CMAKE_SOURCE_DIR}/legacy/src/cycle_over_realizations_undirected.cpp
    ${CMAKE_SOURCE_DIR}/legacy/src/em_update_undirected.cpp
    )

# libraries
add_library(multitensor_default INTERFACE)
target_include_directories(multitensor_default
    INTERFACE
        ${CMAKE_SOURCE_DIR}/legacy/include/
        ${CMAKE_SOURCE_DIR}/legacy/src/)
target_sources(multitensor_default
    INTERFACE
        $<BUILD_INTERFACE:${multitensor_legacy_src} ${multitensor_legacy_default_src}>) # for the IDE to see the library

add_library(multitensor_undirected INTERFACE)
target_include_directories(multitensor_undirected
    INTERFACE
        ${CMAKE_SOURCE_DIR}/legacy/include/
        ${CMAKE_SOURCE_DIR}/legacy/src/)
target_sources(multitensor_undirected
    INTERFACE
        $<BUILD_INTERFACE:${multitensor_legacy_src} ${multitensor_legacy_undirected_src}>) # for the IDE to see the library

add_library(multitensor INTERFACE)
target_include_directories(multitensor
    INTERFACE
        ${CMAKE_SOURCE_DIR}/include/
        ${CMAKE_SOURCE_DIR}/src/)
target_sources(multitensor
    INTERFACE
        $<BUILD_INTERFACE:${multitensor_src} >) # for the IDE to see the library


add_library(multitensor_utils
    ${CMAKE_SOURCE_DIR}/applications/include/app_utils.h
    ${CMAKE_SOURCE_DIR}/applications/src/app_utils.cpp)
target_link_libraries(multitensor_utils PUBLIC multitensor Boost::program_options)
target_include_directories(multitensor_utils
    PUBLIC
    ${CMAKE_SOURCE_DIR}/applications/include/
    ${CMAKE_SOURCE_DIR}/applications/src/)

# executables
add_executable(Multitensor_main legacy/applications/src/MultiTensor_main.cpp)
target_link_libraries(Multitensor_main multitensor_default Boost::system Boost::filesystem Boost::program_options)

add_executable(Multitensor_undirected legacy/applications/src/MultiTensor_undirected.cpp)
target_link_libraries(Multitensor_undirected multitensor_undirected Boost::system Boost::filesystem Boost::program_options)

add_executable(Multitensor applications/src/MultiTensor.cpp)
target_link_libraries(Multitensor multitensor multitensor_utils Boost::system Boost::filesystem Boost::program_options)

set_target_properties(
    Multitensor_main Multitensor_undirected Multitensor
    PROPERTIES
    FOLDER "Applications")

# Unit tests
add_executable(test_multitensor
    tests/test_graph.cpp
    tests/test_multitensor.cpp
    tests/test_solver.cpp
    tests/test_tensor.cpp)
target_link_libraries(test_multitensor multitensor
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
    test_program
    output_check_clusters
    cluster_or_assignments)

    if("${test_program}" STREQUAL "Multitensor_main")
        set(output_var
            -E _main.dat
            )
    else()
        set(output_var
            -E _undirected.dat
        )
    endif()

    add_test(
        NAME "command-line-${test_name}-run"
        COMMAND ${test_program}
            ${ARGN}
            ${output_var}
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/data
            )
endfunction()

add_test_cmd_line(
    "main"
    Multitensor_main
    -k 2
    -l 3
    -a adjacency.dat
    )

add_test_cmd_line(
    "undirected"
    Multitensor_undirected
    -k 2
    -l 3
    -a adjacency.dat
    )


#########################
# Python Bindings + Tests
include(cmake/CMakePython.cmake)
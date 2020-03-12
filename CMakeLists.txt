# Copyright 2019, Max Planck Society.
# Distributed under the "GNU GPL v3" licence.
# (See accompanying file LICENSE.md)

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
    ${CMAKE_SOURCE_DIR}/include/multitensor/parameters.hpp
    ${CMAKE_SOURCE_DIR}/include/multitensor/solver.hpp
    ${CMAKE_SOURCE_DIR}/include/multitensor/tensor.hpp
    ${CMAKE_SOURCE_DIR}/include/multitensor/utils.hpp
    ${CMAKE_SOURCE_DIR}/include/multitensor/main.hpp)
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
    ${CMAKE_SOURCE_DIR}/applications/include/app_utils.hpp
    ${CMAKE_SOURCE_DIR}/applications/src/app_utils.cpp)
target_link_libraries(multitensor_utils PUBLIC multitensor Boost::filesystem)
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
target_link_libraries(Multitensor multitensor multitensor_utils)

set_target_properties(
    Multitensor_main Multitensor_undirected Multitensor
    PROPERTIES
    FOLDER "Applications")

# Unit tests
add_executable(test_multitensor
    tests/fixtures.hpp
    tests/test_graph.cpp
    tests/test_multitensor.cpp
    tests/test_solver.cpp
    tests/test_tensor.cpp
    tests/test_utils.cpp)
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
    test_program)

    # Set test folder and copy adjacency file
    set(test_folder ${CMAKE_BINARY_DIR}/tests/${test_name})
    file(COPY ${CMAKE_SOURCE_DIR}/data/adjacency.dat DESTINATION ${test_folder}/)

    # Run
    add_test(
        NAME "command-line-${test_name}-run"
        COMMAND ${test_program} ${ARGN}
        WORKING_DIRECTORY ${test_folder})

    if("${test_program}" STREQUAL "Multitensor")
        add_test(
            NAME "command-line-${test_name}-verify-content"
            COMMAND bash "-c" "./compare.sh ${test_folder}/results"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/data)

        set_tests_properties(
            "command-line-${test_name}-verify-content"
            PROPERTIES
            DEPENDS "command-line-${test_name}-run")
    endif()

endfunction()

add_test_cmd_line(
    "main"
    Multitensor
    --a adjacency.dat
    --k 2
    --s 5489)

add_test_cmd_line(
    "legacy_main"
    Multitensor_main
    -a adjacency.dat
    -k 2
    -l 3)

add_test_cmd_line(
    "legacy_undirected"
    Multitensor_undirected
    -a adjacency.dat
    -k 2
    -l 3)


#########################
# Python Bindings + Tests
include(cmake/python.cmake)


#########################
# Documentation (Doxygen)
include(cmake/documentation.cmake)
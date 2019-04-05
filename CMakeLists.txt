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
add_definitions(-DMULTI_TENSOR_VERBOSE=0)

#################################
# thirdparties and configurations
include(cmake/compiler_settings.cmake)
include(cmake/thirdparty.cmake)


################################
# Main libraries and executables
set(multitensor_src
    ${CMAKE_SOURCE_DIR}/include/main.hpp
    ${CMAKE_SOURCE_DIR}/include/mlg.hpp
    ${CMAKE_SOURCE_DIR}/src/tools_assortative.cpp
    ${CMAKE_SOURCE_DIR}/src/tools.cpp
    )
set(multitensor_src_default
    ${CMAKE_SOURCE_DIR}/include/cycle_over_realizations.hpp
    ${CMAKE_SOURCE_DIR}/src/cycle_over_realizations.cpp
    ${CMAKE_SOURCE_DIR}/src/em_update.cpp
    )
set(multitensor_src_undirected
    ${CMAKE_SOURCE_DIR}/include/cycle_over_realizations_undirected.hpp
    ${CMAKE_SOURCE_DIR}/src/cycle_over_realizations_undirected.cpp
    ${CMAKE_SOURCE_DIR}/src/em_update_undirected.cpp
    )

add_library(multitensor_default INTERFACE)
target_include_directories(multitensor_default
    INTERFACE
        ${CMAKE_SOURCE_DIR}/include/
        ${CMAKE_SOURCE_DIR}/src/)
target_sources(multitensor_default
    INTERFACE
        $<BUILD_INTERFACE:${multitensor_src} ${multitensor_src_default}>) # for the IDE to see the library

add_library(multitensor_undirected INTERFACE)
target_include_directories(multitensor_undirected
    INTERFACE
        ${CMAKE_SOURCE_DIR}/include/
        ${CMAKE_SOURCE_DIR}/src/)
target_sources(multitensor_undirected
    INTERFACE
        $<BUILD_INTERFACE:${multitensor_src} ${multitensor_src_undirected}>) # for the IDE to see the library

# executables
add_executable(Multitensor_main applications/MultiTensor.cpp)
target_link_libraries(Multitensor_main multitensor_default Boost::system Boost::filesystem Boost::program_options)

add_executable(Multitensor_undirected applications/MultiTensor_undirected.cpp)
target_link_libraries(Multitensor_undirected multitensor_undirected Boost::system Boost::filesystem Boost::program_options)

set_target_properties(
    Multitensor_main Multitensor_undirected
    PROPERTIES
    FOLDER "Applications")

# Unit tests


#########################
# Python Bindings + Tests
include(cmake/CMakePython.cmake)
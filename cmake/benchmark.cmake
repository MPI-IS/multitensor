# Copyright (c) 2019, Max Planck Society / Software Workshop - Max Planck Institute for Intelligent Systems
# Distributed under the GNU GPL license version 3
# See file LICENSE.md or at https://github.com/MPI-IS/multitensor/LICENSE.md

#
# This file contains the logic for running the benchmark
#

if(NOT MULTITENSOR_ENABLE_BENCHMARK)
    return()
endif()

message(STATUS "[BENCHMARK] Configuring benchmark")

# Paths
set(BENCHMARK_SRC_DIR ${CMAKE_SOURCE_DIR}/benchmark)
set(BENCHMARK_TARGET_DIR ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${BENCHMARK_SRC_DIR} DESTINATION ${BENCHMARK_TARGET_DIR})

# Source files
add_custom_target(benchmark_script
    SOURCES
        ${BENCHMARK_SRC_DIR}/benchmark_runner.py)
set_target_properties(benchmark_script
    PROPERTIES
    FOLDER Benchmarks)

# Copyright (c) 2019, Max Planck Society / Software Workshop - Max Planck Institute for Intelligent Systems
# Distributed under the GNU GPL license version 3
# See file LICENSE.md or at https://github.com/MPI-IS/multitensor/LICENSE.md
#
# This file contains the logic for generating the documentation with Doxygen
#

if(NOT DOXYGEN_FOUND)
    return()
endif()

message(STATUS "[DOXYGEN] Configuring the documentation")

# Paths
set(DOXYGEN_SRC_DIR ${CMAKE_SOURCE_DIR}/doc)
set(DOXYGEN_TARGET_DIR ${CMAKE_CURRENT_BINARY_DIR}/doc)

configure_file(
    ${DOXYGEN_SRC_DIR}/Doxyfile.in
    ${DOXYGEN_TARGET_DIR}/Doxyfile)

add_custom_target(doxygen
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE}
    COMMENT "[DOXYGEN] Generating documentation"
    DEPENDS
        ${DOXYGEN_TARGET_DIR}/Doxyfile
        ${multitensor_src}
    WORKING_DIRECTORY
        ${DOXYGEN_TARGET_DIR}
    SOURCES
        ${CMAKE_SOURCE_DIR}/cmake/doxygen_doc.cmake
        ${DOXYGEN_SRC_DIR}/Doxyfile.in
        ${CMAKE_SOURCE_DIR}/README.md
)

set_target_properties(doxygen
    PROPERTIES
        FOLDER "Documentation")

#
# This file contains the logic for generating the documentation with Doxygen
#

if(NOT DOXYGEN_FOUND)
    return()
endif()

message(STATUS "[DOC] Configuring the documentation")

# Paths
set(DOXYGEN_SRC_DIR ${CMAKE_SOURCE_DIR}/doc)
set(DOXYGEN_TARGET_DIR ${CMAKE_CURRENT_BINARY_DIR}/doc)

configure_file(
    ${DOXYGEN_SRC_DIR}/Doxyfile.in
    ${DOXYGEN_TARGET_DIR}/Doxyfile)

add_custom_target(Doxygen
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE}
    COMMENT "Generating Doxygen documentation"
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

set_target_properties(Doxygen
    PROPERTIES
        FOLDER "Documentation")

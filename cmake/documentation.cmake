#
# This file contains the logic for generating the documentation with Doxygen
#

if(NOT DOXYGEN_FOUND)
    return()
endif()

message(STATUS "[DOC] Configuring the documentation")

configure_file(
    ${CMAKE_SOURCE_DIR}/doc/Doxyfile.in
    ${CMAKE_CURRENT_BINARY_DIR}/documentation/Doxyfile)

add_custom_target(Doxygen
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE}
    COMMENT "Generating Doxygen documentation"
    DEPENDS
        ${CMAKE_CURRENT_BINARY_DIR}/documentation/Doxyfile
        ${multitensor_src}
    WORKING_DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR}/documentation
    SOURCES
        ${CMAKE_SOURCE_DIR}/cmake/documentation.cmake
        ${CMAKE_SOURCE_DIR}/doc/Doxyfile.in
        ${CMAKE_SOURCE_DIR}/README.md
)

set_target_properties(Doxygen
    PROPERTIES
        FOLDER "Documentation")

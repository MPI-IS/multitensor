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
set(_cython_generated_file "${CMAKE_CURRENT_BINARY_DIR}/multitensor.cpp")
if(CYTHON_PROGRAM)
    add_custom_command(
        OUTPUT
            ${_cython_generated_file}
        COMMAND
            ${CYTHON_PROGRAM}
        ARGS
            --cplus -3
            "${PYTHON_PACKAGE_SRC_DIR}/multitensor.pyx"
            -o
            ${_cython_generated_file}
        COMMENT
            "[PYTHON] Generating the binding files with Cython"
        DEPENDS
            "${PYTHON_PACKAGE_SRC_DIR}/multitensor.pyx"
            ${multitensor_src})
else()
    message(STATUS "[Python] Copy the .ccp file form the repository")
    configure_file(
        ${PYTHON_PACKAGE_SRC_DIR}/multitensor.cpp
        ${_cython_generated_file}
        COPYONLY)
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
add_test(
    NAME multitensor_python_test
    COMMAND ${Python3_EXECUTABLE} -m unittest discover -vvv
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/python)
# Add path to the library
set_tests_properties(multitensor_python_test
    PROPERTIES
        ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:multitensor_py>")

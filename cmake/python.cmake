###########################
# Python Bindings + Tests
###########################

if(MULTI_TENSOR_PYTHON_EXTENSIONS)

    message(STATUS "[PYTHON] Configuring Python bindings")

    # Paths
    set(MULTI_TENSOR_PYTHON_ROOT ${CMAKE_SOURCE_DIR}/python)

    # Gather python files
    file(GLOB PythonFiles ${MULTI_TENSOR_PYTHON_ROOT}/*.py)

    add_custom_target(multi_tensor_py ALL)
    foreach(python_file ${PythonFiles})
        add_custom_command(
            TARGET multi_tensor_py
            PRE_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy ${python_file} ${CMAKE_CURRENT_BINARY_DIR})
    endforeach()

    # Functional tests
    add_test(
        NAME "command-line-python-run"
        COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_BINARY_DIR}/main.py"
            -a=main/adjacency.dat
            -k=2
            -l=4
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/data
            )

    #[[
    # Sphinx Documentation for Python extension
    file(COPY ${ACC_KMEANS_PYTHON_ROOT}/python_doc DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    add_custom_target(Sphinx
        COMMAND make html
        COMMENT "Generating Sphinx documentation"
        DEPENDS ${mpi_kmeans_py}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/python_doc
        )
    add_dependencies(Sphinx mpi_kmeans_py)
    ]]
endif()

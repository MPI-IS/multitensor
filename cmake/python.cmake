###########################
# Python Bindings + Tests
###########################

if(MULTI_TENSOR_PYTHON_EXTENSIONS)
    message(STATUS "[PYTHON] Configuring Python bindings")

    # Paths
    set(MULTI_TENSOR_PYTHON_ROOT ${CMAKE_SOURCE_DIR}/python)

    # We need the boost library path for building the wheel
    set(BOOST_LIB_DIR $<$<CONFIG:Debug>:${Boost_LIBRARY_DIR_DEBUG}>$<$<CONFIG:RelWithDebInfo>:${Boost_LIBRARY_DIR_RELEASE}>$<$<CONFIG:Release>:${Boost_LIBRARY_DIR_RELEASE}>$<$<CONFIG:MinSizeRel>:${Boost_LIBRARY_DIR_RELEASE}>)

    # Build a wheel
    add_custom_target(multitensor_py ALL)
    add_custom_command(
      TARGET multitensor_py
      POST_BUILD
      WORKING_DIRECTORY ${MULTI_TENSOR_PYTHON_ROOT}
      COMMAND ${Python3_EXECUTABLE} -m pip install wheel
      COMMAND ${Python3_EXECUTABLE} -m pip install cython
      COMMAND CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}
              BOOST_INCLUDE_DIR=${Boost_INCLUDE_DIR} BOOST_LIBRARY_DIR=${BOOST_LIB_DIR}
              ${Python3_EXECUTABLE} setup.py bdist_wheel)
    add_custom_command(
      TARGET multitensor_py
      POST_BUILD
      WORKING_DIRECTORY ${MULTI_TENSOR_PYTHON_ROOT}
      COMMAND ${Python3_EXECUTABLE} setup.py clean --all)

    # Copy the wheel into a binary directory
    file(GLOB python_wheel ${MULTI_TENSOR_PYTHON_ROOT}/dist/*.whl)
    file(COPY ${python_wheel} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

    # This is a target that simply installs nose and the package
    # itself just before the tests.
    add_custom_target(multitensor_py_unittest_prepare ALL)
    add_custom_command(
      TARGET multitensor_py_unittest_prepare
      POST_BUILD
      WORKING_DIRECTORY ${MULTI_TENSOR_PYTHON_ROOT}
      COMMAND ${Python3_EXECUTABLE} -m pip install nose)
    add_custom_command(
      TARGET multitensor_py_unittest_prepare
      POST_BUILD
      WORKING_DIRECTORY ${MULTI_TENSOR_PYTHON_ROOT}
      COMMAND CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}
              BOOST_INCLUDE_DIR=${Boost_INCLUDE_DIR} BOOST_LIBRARY_DIR=${BOOST_LIB_DIR}
              ${Python3_EXECUTABLE} setup.py install)

    #[[
    # Sphinx Documentation for Python extension
    file(COPY ${ACC_KMEANS_PYTHON_ROOT}/python_doc DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    add_custom_target(Sphinx
        COMMAND make html
        COMMENT "Generating Sphinx documentation"
        DEPENDS ${mpi_kmeans_py}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/python_doc)
    add_dependencies(Sphinx mpi_kmeans_py)
    ]]

    # Functional tests.
    add_test(
      NAME multitensor_py_unittest
      # We work from the outside directory so that `multitensor` (the
      # original source directory) is not in `PYTHONPATH`. Otherwise it
      # shadows the package installed in virtual environment.
      COMMAND ${Python3_EXECUTABLE} -m nose --no-path-adjustment --with-xunit --verbose ${MULTI_TENSOR_PYTHON_ROOT}/tests.py)
endif()

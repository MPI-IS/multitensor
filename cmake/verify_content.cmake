# Cmake script to verify the content of a run. Variables are
# - 'results' points to a folder containing the results of the tests
# - 'ground_truth' points to a folder containing the ground truth

set(DATALIST u v w)
foreach(X IN LISTS DATALIST)
    FILE(GLOB file1 ${results}/${X}*.dat)
    FILE(GLOB file2 ${ground_truth}/${X}*_compare.txt)
    if(NOT "${file2}" STREQUAL "") # allows to skip v for undirected case
        execute_process(
            COMMAND ${CMAKE_COMMAND} -E compare_files ${file1} ${file2}
            RESULT_VARIABLE _not_the_same
        )
        if(${_not_the_same})
            message(FATAL_ERROR "Verification for ${X} failed.")
        else()
            message(STATUS "Verification for ${X} passed.")
        endif()
    endif()
endforeach()

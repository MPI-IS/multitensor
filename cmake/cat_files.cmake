#
# This files helps to concatenate two files
#

file(READ "${file1}" FILE1)
file(WRITE "${output}" "${FILE1}")
file(READ "${file2}" FILE2)
file(APPEND "${output}" "${FILE2}")

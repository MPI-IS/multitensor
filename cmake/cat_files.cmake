# Copyright (c) 2019, Max Planck Society / Software Workshop - Max Planck Institute for Intelligent Systems
# Distributed under the GNU GPL license version 3
# See file LICENSE.md or at https://github.com/MPI-IS/multitensor/LICENSE.md
#
# This files helps to concatenate two files
#

file(READ "${file1}" FILE1)
file(WRITE "${output}" "${FILE1}")
file(READ "${file2}" FILE2)
file(APPEND "${output}" "${FILE2}")

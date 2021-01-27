#!/bin/bash

# Copyright (c) 2019, Max Planck Society / Software Workshop - Max Planck Institute for Intelligent Systems
# Distributed under the GNU GPL license version 3
# See file LICENSE.md or at https://github.com/MPI-IS/multitensor/LICENSE.md

# Small script comparing the different output files produced by the algorithm

for i in u v w;
do
    if [ -e ${i}_*_compare.txt ]; then
        diff ${i}_*_compare.txt $1/${i}_out.dat
    fi
done

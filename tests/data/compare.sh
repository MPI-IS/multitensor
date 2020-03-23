#!/bin/bash

# Scmall scripts comparing the different output files produced by the algorithm
for i in u v w;
do
    if [ -e ${i}_*_compare.txt ]; then
        diff ${i}_*_compare.txt $1/${i}_out.dat
    fi
done
#!/bin/bash

# Scmall scripts comparing the different output files produced by the algorithm
for i in u v w;
do
    diff ${i}_out.dat $1/${i}_K2_compare.txt;
done
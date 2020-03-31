#!/bin/bash

# Scmall scripts comparing the different output files produced by the algorithm
for i in u v w;
do
    diff ${i}_K2_compare.txt $1/${i}_out.dat ;
done
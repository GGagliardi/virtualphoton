#!/bin/bash
export OMP_SET_NUM_THREADS=5
g++ -o virt_mod.o virtual_model.cpp -g  -Wall -std=c++17 -O3 -rdynamic -L. -lmpfr -fopenmp

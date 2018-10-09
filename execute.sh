#!/bin/bash -x

export LD_LIBRARY_PATH="$PWD/lib:$LD_LIBRARY_PATH"

mpirun -np 6 ./montecarlo_x86_64_c++ 

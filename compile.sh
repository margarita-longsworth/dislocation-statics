#!/bin/bash -x

## GCC ParaStationMPI toolchain

#module load GCC
#module load ParaStationMPI
#module load CMake

## Intel IntelMPI toolchain

#module load Intel
#module load IntelMPI
#module load CMake

cd 0K_src_a2951c0
make clean
cd ..

mkdir build
cd build
cmake ../monte_carlo 
make
cp montecarlo_x86_64_c++ ..
cd ..


#!/bin/bash -x
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=00:30:00
#SBATCH --partition=devel

export LD_LIBRARY_PATH="$PWD/lib:$LD_LIBRARY_PATH"

## GCC ParaStationMPI toolchain

module load GCC
module load ParaStationMPI

srun ./bin/montecarlo_x86_64_c++

## Intel IntelMPI toolchain

#module load Intel
#module load IntelMPI

#srun ./bin/montecarlo_x86_64_c++


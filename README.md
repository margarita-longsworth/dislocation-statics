# Monte Carlo MD Simulation

This repository contains the Monte Carlo Simulation of Margartia Longsworth.

# Setup

The simulation code uses the CMake buildsystem generator to setup up a
Make-centric build system.

## Prerequesites

The simulation uses the [ITAP MD simulation library](https://github.com/itapmd/imd).
The CMake setup distinguishes between the GNU and Intel compiler family and
will internally trigger the build of the ITAP MD library with the corresponding
compiler family. However, as the ITAP MD library is not integrated into the
CMake build, an in-source build is performed. Therefore, you need to ensure the
library directory is clean, before you initially setup and build the Monte
Carlo simulation code.

## Creating the Makefiles

It is recommended to perform an out-of-source build, to separate build files
from source files. For this you can create an empty build directory in the
project root.

```sh
% mkdir _build
% cd _build/
```

You can then generate the Makefiles necessary to build the code using the
`cmake` command and choosing the appropriate toolchain out of
`cmake/toolchains`.

Below, you find some example for available toolchains:

* Intel compiler & Intel MPI

```sh
% cmake ../monte_carlo \
        -DCMAKE_TOOLCHAIN_FILE=Intel-IntelMPI.cmake
```
* Intel compiler & ParaStation MPI

```sh
% cmake ../monte_carlo \
        -DCMAKE_TOOLCHAIN_FILE=Intel-ParaStation.cmake
```

* GNU compiler & ParaStation MPI

```sh
% cmake ../monte_carlo \
        -DCMAKE_TOOLCHAIN_FILE=GNU-ParaStation.cmake
```

After the build system is set up, you can execute the make process.

```sh
% make
```

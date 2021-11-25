# SPM: Finite Volume Solver for Radiative Transfer

[![CI](https://github.com/CSMMLab/KiT-RT/actions/workflows/c-cpp.yml/badge.svg)](https://github.com/CSMMLab/KiT-RT/actions/workflows/c-cpp.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

The SPM framework is a high-performance open source platform for solving the linear Boltzmann equation, which is commonly adopted to describe radiative transfer.
The implementation solvers includes SN (discrete ordinate), PN (spherical harmonics), and MN (moment closure) methods.
A short description of kinetic theory can be found [here](https://kit-rt.readthedocs.io/en/develop/physics.html).

## Build
### Required dependencies
 - Compiler with C++17 support
 - cmake >= v3.16
 - LAPACK
 - OpenMP
 - MPI
 - python3
 - VTK
 - git

### Obtain submodules
Note that an **active internet connection is required for the first build** in order to download the suitable versions of the required submodules!
For the first build only, download all submodules:

```bash
git submodule update --init --recursive
```

### Compile the code
The executable file can be build with cmake workflow:
 
```bash 
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
make -j N
```

## Run
Execute the compiled binary by handing over a [valid config file](https://kit-rt.readthedocs.io/en/latest/configFiles.html), e.g.:

```bash
./SPM ../input/linesource_SN.cfg
```

In order to run the code in parallel execute:

```bash
OMP_NUM_THREADS=N mpirun -np J ./SPM ../examples/linesource_SN.cfg
```

with `N` equal to the number of shared memory threads and `J` equal to the number of distrubuted memory threads.

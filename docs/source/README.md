 [![Build Status](https://travis-ci.org/DavidAce/DMRG.svg?branch=master)](https://travis-ci.org/DavidAce/DMRG)



# DMRG++

## Introduction
  [Density matrix renormalization group](https://en.wikipedia.org/wiki/Density_matrix_renormalization_group) (DMRG) is a variational numerical technique used to study many-body
  quantum systems. It works by optimizing a trial wave function in the form of a [Matrix Product States](https://en.wikipedia.org/wiki/Matrix_product_states) (MPS), to find either the
  groundstate or an eigenstate of a 1D quantum system with high precision. DMRG++ includes 4 different MPS-based algorithms for 1D systems:

  - ***i*DMRG:** *infinite* DMRG. Finds the groundstate of infinite and translationally invariant systems.
  - ***f*DMRG:** *finite* DMRG. Finds the groundstate of finite systems, not necessarily translationally invariant.
  - ***x*DMRG:** *Excited state* DMRG. Finds highly excited (mid-spectrum) eigenstates of finite systems.
  - ***i*TEBD:** *Imaginary Time Evolving Block Decimation*. Finds the ground state of infinite and translationally invariant systems using unitary operators that perform imaginary time evolution.

The program is controlled through an input file, whose path (full or relative to the binary) is given as input argument in the command line. See the [Installation](#installation) section below.

Included here are two 1D models of spin chains, the *Quantum Ising model with transverse field* and the *Self-dual quantum Ising model with random couplings and random fields*. The choice of model
is done in the input configuration file.

### Working Notes (in construction)

 Go to the [working notes](https://github.com/DavidAce/Notebooks/blob/master/DMRG%2B%2B/DMRG%2B%2B.pdf) on the theoretical aspects of this implementation.



---
## Installation
### Quick start
- Git clone and build with `./build.sh`.
- Run with `./run.sh` after setting your options in `./input/input.cfg`.
- Profit from `./output/output.h5`.

### Build
Git clone or copy & extract the project into a folder of your choosing.
**Make sure there are no spaces in your path!**.

The project can be built with a single command from a unix terminal.
Simply launch the script `.\build.sh` found in the root folder to trigger a CMake build.

The script takes optional arguments, run `.\build.sh -h` to learn more.

**Alternatively**, if you intend to develop or study the source code, some IDE's with CMake support can self-configure from the file CMakeLists.txt found in the project root folder. This
is perhaps an even simpler approach. Recommended: [CLion](https://www.jetbrains.com/clion/download) or [Visual Studio Code](https://code.visualstudio.com/) with the C++ and CMake Tools extensions.


CMake will check for dependencies in the host system. If not found, it will download and install these automatically to a folder `libs` in the project root (see *Optional Requirements* below).
If the dependencies are successfully found or installed, the project is built and an executable is generated in `./build/Release/DMRG++`.


### Execution
To run the executable, launch `.\run.sh`, or `.\run.sh -h` to learn more.
If no configuration file is given as argument to `run.sh`, the executable will look for a configuration file located in `./input/input.cfg`.

#### Configuration file
The default configuration file in `./input/input.cfg` contains run-time instructions for DMRG++. Here you can choose the type of simulation, type of model, model parameters,
system size, precision as well as settings for profiling, data storage and console verbosity. Read the comments in the file to learn more.


#### Output file
After execution the results are stored a binary file in HDF5 format. Its location is specified in the configuration file `./input/input.cfg`.
By default this should be in `.output/output.h5`. This file will contain values like the final energies, entanglement entropies, entanglement spectrums and
optionally the final state in MPS form.

To view the data you can use any hdf5-viewer, such as HDFCompass.

The script `analysis/data_analysis.py` (in progress) shows how to analyze the simulation data stored in the hdf5 files. You need to install the python package
`h5py` from pip or conda to read files in the HDF5 format.


### Minimum Requirements
The following software is required to build the project:
 - C++ compiler with support for C++17, OpenMP and libstdc++ standard library implementation  (version >= 7). Tested with two compilers:
    - GNU GCC versions 7 and 8 (these bundle libstdc++)
    - Clang version >= 7.0. (you need to manually install libstdc++ version >= 7, that comes bundled with gcc, for instance from `ppa:ubuntu-toolchain-r/test`)
 - CMake version >= 3.9. If you compile CMake from source, remember to enable `curl` (`./bootstrap --system-curl`). 
 - *(For library compilation only) Fortran compiler, tested with gfortran-7 and gfortran-8.*

**Ubuntu** 17 or higher will have the above versions in the default repositories. For older distributions, use the ppa `ubuntu-toolchain-r/test` to get newer versions.

**Mac OSX** users are advised to use GNU GCC version 7 or 8 from homebrew. Install with `brew install gcc`. Clang from llvm 6.0 might work but you will have to link to GNU's `libstdc++.so` or `libstdc++.a` manually. The AppleClang compiler is not supported at all. 


### Optional Requirements
The compilation of DMRG++ requires several libraries. To meet the requirements, you have two options:

    1. **Automatic**: let CMake download and compile the libraries below from source into a local folder `./libs`. This is the default behavior if the library is not found on your system. Note that this does *NOT* make a system-wide install, so you can safely delete the `./libs` folder.
    2. **Manual**: install the libraries yourself with your favorite package manager (e.g. `apt` in ubuntu or `brew` in OSX). The build script will always attempt to find the libraries in your system first.
    3. **Manual with environment modules** (in construction). You can also load *modules* from the ubuntu command-line tool *environment-modules*. Simply export the variable `<LIBNAME>_DIR` to let CMake know where to look.

 If the compilation halts due to any library failing to compile or link, you can try installing/uninstalling that library from your package manager.

#### List of libraries

 - **BLAS** and **LAPACK**. Required for Arpack. You can choose either [Intel MKL](https://software.intel.com/en-us/mkl) or [OpenBLAS](https://github.com/xianyi/OpenBLAS). If not found, OpenBLAS is downloaded automatically. Note that OpenBLAS requires Fortran to compile from source. If both MKL and OpenBLAS are found in the system, MKL is preferred.
 - [**Eigen**](http://eigen.tuxfamily.org) for tensor and matrix and linear algebra (tested with version >= 3.3).
 - [**Arpack**](https://github.com/opencollab/arpack-ng) Eigenvalue solver based on Fortran. Note that this in turn requires LAPACK and BLAS libraries, both of which are included in OpenBLAS.
 - [**Arpackpp**](https://github.com/m-reuter/eigsolver_properties) C++ frontend for Arpack.
 - [**h5pp**](https://github.com/DavidAce/h5pp) a wrapper for HDF5.
 - [**HDF5**](https://support.hdfgroup.org/HDF5/) for hdf5 binary output file support (tested with version >= 1.10).
 - [**spdlog**](https://github.com/gabime/spdlog) logger).

---


## License
Open source under [MPL2](https://www.mozilla.org/MPL/2.0/).

## Author
David Aceituno
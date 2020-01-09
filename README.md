[![Build Status](https://travis-ci.org/DavidAce/DMRG.svg?branch=master)](https://travis-ci.org/DavidAce/DMRG)
[![Build Status](https://github.com/DavidAce/DMRG/workflows/Actions/badge.svg)](https://github.com/DavidAce/DMRG/actions)


# DMRG++

## Introduction
  [Density matrix renormalization group](https://en.wikipedia.org/wiki/Density_matrix_renormalization_group) (DMRG) is a variational technique used to study 1D quantum systems. It works by iteratively optimizing a trial wave function in the form of a [Matrix Product State](https://en.wikipedia.org/wiki/Matrix_product_states) (MPS), until obtaining an eigenstate of the system with high precision. DMRG++ includes 4 different algorithms:

  - ***i*DMRG:** *infinite* DMRG. Finds the groundstate of infinite and translationally invariant systems.
  - ***f*DMRG:** *finite* DMRG. Finds the groundstate of finite systems, not necessarily translationally invariant.
  - ***x*DMRG:** *Excited state* DMRG. Finds highly excited (mid-spectrum) eigenstates of finite systems.
  - ***i*TEBD:** *Imaginary Time Evolving Block Decimation*. Finds the ground state of infinite and translationally invariant systems using unitary operators that perform imaginary time evolution.


### Documentation (in construction)
 Go to the [Documentation](https://davidace.github.io/DMRG/) page generated by [Doxygen](http://www.stack.nl/~dimitri/doxygen/).

### Working Notes (in construction)
 Go to the [working notes](https://github.com/DavidAce/Notebooks/blob/master/DMRG%2B%2B/DMRG%2B%2B.pdf) on the theoretical aspects of this implementation.



---
## Installation
### Quick start
- Git clone and build with `./build.sh` (see options with `./build.sh --help`).
- Run with `./build/Release/DMRG++ -i input/input.cfg` after setting your options in `input/input.cfg`.
- Profit from `output/output.h5`.

### Build
Git clone or copy & extract the project into a folder of your choosing.
**Make sure there are no spaces in your path!**.

The project can be built with a single command from a unix terminal.
Simply launch the script `build.sh` found in the root folder to trigger a CMake build.

The script takes optional arguments, run `build.sh --help` to learn more.

CMake will check for dependencies in the host system. If not found, it will download and install these automatically to a folder `libs` in the project root (see *Optional Requirements* below).
If the dependencies are successfully found or installed, the project is built and an executable is generated in `build/Release/DMRG++`.

### Usage
The input configuration file `input/input.cfg` controls the properties of the simulation. `DMRG++` admits custom input files from command-line arguments, e.g. `./build/Release/DMRG++ -i ../CustomFolder/custom.cfg`.

If no configuration file is given as argument, the default is to look for a configuration file located in `input/input.cfg` relative to the project root folder.

#### Configuration file
The default configuration file in `input/input.cfg` contains run-time instructions for DMRG++. Here you can choose the type of simulation, type of model, model parameters,
system size, precision as well as settings for profiling, data storage and console verbosity. Read the comments in the file to learn more.

#### Models
Two models of quantum systems are included, the *(Random) Transverse-field Ising model* and the *Self-dual (random) transverse-field Ising model* located under `source/model`. To add another model, one currently
needs to hand-craft the *Matrix Product Operator* (MPO) of the model, and add it to the source code via a derived class just as the preexisting ones. 
Once implemented, the model is selected using the input configuration file in `input/input.cfg`, with the option `model::model_type `.

#### Output file
After execution the results are stored a binary file in HDF5 format. Its location is specified in the configuration file `input/input.cfg`.
By default this should be in `output/output.h5`. This file will contain values like the final energies, entanglement entropies, entanglement spectrums and
optionally the final state in MPS form.

To view the data you can use any hdf5-viewer, such as HDFCompass.

The script `analysis/data_analysis.py` (in progress) shows how to analyze the simulation data stored in the output files. You need to install the python package
`h5py` from pip or conda to read files in the HDF5 format.


### Minimum Requirements
The following software is required to build the project:
 - C++ compiler with support for C++17 and libstdc++ standard library implementation  (version >= 7). Tested with two compilers:
    - GNU GCC versions 7 and 8 (these bundle libstdc++)
    - Clang version >= 6.0. (you need to manually install libstdc++ version >= 7, that comes bundled with gcc, for instance from `ppa:ubuntu-toolchain-r/test`)
 - CMake version >= 3.10. 
 - *(For library compilation only) Fortran compiler, tested with gfortran-7 and gfortran-8.*
 
**Ubuntu** 17 or higher will have the above versions in the default repositories. For older distributions, use the ppa `ubuntu-toolchain-r/test` to get newer versions.

**Mac OSX** users are advised to use GNU GCC version 7 or 8 from homebrew. Install with `brew install gcc`. Clang from llvm 6.0 might work but you will have to link to GNU's `libstdc++.so` or `libstdc++.a` manually. The AppleClang compiler is not supported at all. 


### Optional Requirements
The compilation of DMRG++ requires several libraries. To meet the requirements, you have two options:

  1. **Automatic**: let CMake download and compile the libraries below from source into a local folder `libs-release` or `libs-debug`. This is an opt-in behavior if the library is not found on your system. Note that this does *NOT* make a system-wide install, so you can safely delete the `libs-<config>` folders.
  2. **Manually**: install the libraries yourself with your favorite package manager (e.g. `conda`,`apt` in ubuntu or `brew` in OSX). The build attempts to find libraries in your local system. 
  3. **Manual with modules from [Easybuild](https://easybuild.readthedocs.io/en/latest/)** (in construction). You can also load *modules* from the ubuntu command-line tool *environment-modules* or *Lmod*.  CMake will look for environment variables such as `EBROOT<libname>` that are defined when loading Easybuild modules.
 
 If the compilation halts due to any library failing to compile or link, you can try installing/uninstalling that library from other sources package manager, or select conan as the preferred download method. This is done
 with the flag `./build.sh --download-method=conan` or directly as a CMake cli parameter `-DDOWNLOAD_METHOD=conan`. This requires conan to be installed in your system e.g., through apt/pip/conda.
 
#### List of libraries
 
 - **BLAS** and **LAPACK**. Required for Arpack. You can choose either [Intel MKL](https://software.intel.com/en-us/mkl) or [OpenBLAS](https://github.com/xianyi/OpenBLAS). If not found, OpenBLAS is downloaded automatically. Note that OpenBLAS requires Fortran to compile from source. If both MKL and OpenBLAS are found in the system, MKL is preferred.
 - [**Eigen**](http://eigen.tuxfamily.org) for tensor and matrix and linear algebra (tested with version >= 3.3).
 - [**Arpack**](https://github.com/opencollab/arpack-ng) Eigenvalue solver based on Fortran. Note that this in turn requires LAPACK and BLAS libraries, both of which are included in OpenBLAS.
 - [**Arpackpp**](https://github.com/m-reuter/eigsolver_properties) C++ frontend for Arpack.
 - [**h5pp**](https://github.com/DavidAce/h5pp) a wrapper for HDF5. 
 - [**spdlog**](https://github.com/gabime/spdlog) A fast logger based on fmt.
 - [**HDF5**](https://support.hdfgroup.org/HDF5/) for output binary output file support (tested with version >= 1.10).
 - [**ceres**](http://ceres-solver.org/) Optimization library. Here we use the L-BFGS routines. 

---

 
## License
Open source under [MPL2](https://www.mozilla.org/MPL/2.0/).

## Author
David Aceituno
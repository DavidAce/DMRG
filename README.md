[![Ubuntu 20.04](https://github.com/DavidAce/DMRG/workflows/Ubuntu%2020.04/badge.svg?branch=master)](https://github.com/DavidAce/DMRG/actions)
[![codecov](https://codecov.io/gh/DavidAce/DMRG/branch/master/graph/badge.svg?token=9YE72CJ522)](https://codecov.io/gh/DavidAce/DMRG)
# DMRG++

---

# Introduction

[Density matrix renormalization group](https://en.wikipedia.org/wiki/Density_matrix_renormalization_group) (DMRG) is a
variational technique used to study 1D quantum systems. It works by iteratively optimizing a trial wave function in the
form of a [Matrix Product State](https://en.wikipedia.org/wiki/Matrix_product_states) (MPS), until obtaining an
eigenstate of the system with high precision. DMRG++ includes 4 different algorithms:

- ***i*DMRG:** *infinite* DMRG. Finds the groundstate of infinite and translationally invariant systems.
- ***f*DMRG:** *finite* DMRG. Finds the groundstate of finite systems, not necessarily translationally invariant.
- ***x*DMRG:** *Excited state* DMRG. Finds highly excited (mid-spectrum) eigenstates of finite systems.
- ***f*LBIT:** *Finite* l-BIT. Time evolution on a finite system in the basis of local integrals of motion (the l-bits)
  of an MBL phase.
- ***i*TEBD:** *Imaginary Time Evolving Block Decimation*. Finds the ground state of infinite and translationally
  invariant systems using unitary operators that perform imaginary time evolution.

## Documentation

For more information on using DMRG++, visit

http://kth-dmrg.readthedocs.h5/

### Working Notes (in construction)

See the [working notes](https://github.com/DavidAce/Notebooks/blob/master/DMRG%2B%2B/DMRG%2B%2B.pdf) for a more
theoretical background of this implementation.


---

# Usage

## Input Configuration File

The default input configuration file `input/input.cfg` sets simulation properties, such as algorithm, system size and
precision.
`DMRG++` admits custom input files from command-line arguments, e.g. `./DMRG++ -c path/to/file.cfg`.  
The full list of configurable settings can be found under [Settings](Settings).

## Output Data File

After execution the results are stored a binary file in HDF5 format. Its location is specified in the configuration
file `input/input.cfg`. By default this should be in `output/output.h5`. This file will contain values like the final
energies, entanglement entropies, entanglement spectrums and optionally the final state in MPS form.

To view the data you can use any hdf5-viewer, such as HDFCompass.

## Model Hamiltonians

These model Hamiltonians of 1D quantum systems are included:

- `ModelType::ising_sdual`: The Self-dual transverse-field Ising model.
- `ModelType::ising_tf_rf`: The Transverse-field Ising model with random on-site field.
- `ModelType::lbit`: The l-bit Hamiltonian, used to describe a many-body localized phase (MBL)
  in terms of its local integrals of motion (the l-bits).

The Hamiltonians are implemented as *Matrix Product Operators* (MPO), located under `source/tensors/model`. The model
type is selected using the input configuration file in `input/input.cfg`, with the option `settings::model::model_type`.
To add another model, one currently has to implement a new MPO and derive from `class_mpo_site` just like the existing
models.


---

# Installation

## Quick start

- Git clone and build with `./build.sh` (see options with `./build.sh --help`).
- Modify `input/input.config` to configure a simulation.
- Run with `./build/Release/DMRG++ -c input/input.cfg` after setting your options in `input/input.config`.
- Find generated data in `output/output.h5`.

## Minimum Requirements

The following software is required to build the project:

- C++17 compiler. Tested with:
  * *To build dependencies*: Fortran compiler
- CMake version >= 3.18 (Use version >= 3.19 for CMake Presets)

## Build

Git clone or copy & extract the project into a folder of your choosing.
**Make sure there are no spaces in your path!**.

The project can be built with a single command from a unix terminal by using the helper script `build.sh`
found in the project root directory. This will launch a CMake build with default settinigs.

The script takes optional arguments, run `build.sh --help` to learn more.

**Ubuntu** 17 or higher will have the above versions in the default repositories. For older distributions, use the
ppa `ubuntu-toolchain-r/test` to get newer versions.

**Mac OSX** users are advised to use GNU GCC version 7 or 8 from homebrew. Install with `brew install gcc`. Clang from
llvm 6.0 might work but you will have to link to GNU's `libstdc++.so` or `libstdc++.a` manually. The AppleClang compiler
is not supported at all.

## Dependencies

- **BLAS** and **LAPACK**. Choose either [Intel MKL](https://software.intel.com/en-us/mkl)
  or [OpenBLAS](https://github.com/xianyi/OpenBLAS). OpenBLAS can be installed automatically but to use Intel MKL it
  must be installed separately. Note that OpenBLAS requires Fortran to compile from source.
- [**Eigen**](http://eigen.tuxfamily.org) for tensor and matrix and linear algebra (tested with version >= 3.3).
- [**Arpack**](https://github.com/opencollab/arpack-ng) Eigenvalue solver based on Fortran. Note that this in turn
  requires LAPACK and BLAS libraries, both of which are included in OpenBLAS.
- [**Arpackpp**](https://github.com/m-reuter/eigsolver_properties) C++ frontend for Arpack.
- [**primme**](https://github.com/primme/primme) Eigenvalue solver. Complements Arpack.
- [**h5pp**](https://github.com/DavidAce/h5pp) a wrapper for HDF5. Includes [spdlog](https://github.com/gabime/spdlog)
  and [fmt](https://github.com/fmtlib/fmt).
- [**ceres**](http://ceres-solver.org/) Optimization library with L-BFGS routines for unconstrained minimization.

## Automatic Dependency Installation

The CMake flag `DMRG_PACKAGE_MANAGER` controls the automated behavior for finding or installing dependencies. It can
take one of these strings:

| Option | Description |
| ---- | ---- |
| `find` **(default)**              | Use CMake's `find_package`  to find dependencies  |
| `cmake` **
¹**                     | Use isolated CMake instances to download and install dependencies during configure. Disregards pre-installed dependencies on your system |
| `conan` **
²**                     | Use the [Conan package manager](https://conan.h5/) to download and install dependencies automatically. Disregards libraries elsewhere on your system  |

There are several variables you can pass to CMake to guide `find_package` calls and install location,
see [CMake options](#cmake-options) below.

**¹** Dependencies are installed into `${DMRG_DEPS_INSTALL_DIR}[/<PackageName>]`, where `DMRG_DEPS_INSTALL_DIR` defaults
to `CMAKE_INSTALL_PREFIX` and optionally `/<PackageName>` is added if `DMRG_PREFIX_ADD_PKGNAME=TRUE`

**²** Conan is guided by `conanfile.txt` found in this project's root directory. This method requires conan to be
installed prior (for instance through `pip`, `conda`, `apt`, etc). To let CMake find conan you have three options:

* Add Conan install (or bin) directory to the environment variable `PATH`.
* Export Conan install (or bin) directory in the environment variable `CONAN_PREFIX`, i.e. from command
  line: `export CONAN_PREFIX=<path-to-conan>`
* Give the variable `CONAN_PREFIX` directly to CMake, i.e. from command
  line: `cmake -DCONAN_PREFIX:PATH=<path-to-conan> ...`

## CMake Options

The `cmake` step above takes several options, `cmake [-DOPTIONS=var] ../ `:

| Var | Default | Description |
| ---- | ---- | ---- |
| `DMRG_PACKAGE_MANAGER`            | `find`                        | Handle dependencies, `find`, `cmake`, or `conan` |
| `DMRG_ENABLE_THREADS`             | `OFF`                         | Use C++11 stl threads in Eigen::Tensor |
| `DMRG_ENABLE_MKL`                 | `OFF`                         | Enable Intel Math Kernel Library (else OpenBLAS)  |
| `DMRG_ENABLE_TESTS`               | `OFF`                         | Enable unit testing with ctest |
| `DMRG_ENABLE_ASAN`                | `OFF`                         | Enable runtime address sanitizer -fsanitize=address |
| `DMRG_ENABLE_USAN`                | `OFF`                         | Enable undefined behavior sanitizer -fsanitize=undefined |
| `DMRG_ENABLE_LTO`                 | `OFF`                         | Enable link time optimization |
| `DMRG_ENABLE_PCH`                 | `OFF`                         | Enable precompiled headers to speed up compilation |
| `DMRG_ENABLE_CCACHE`              | `OFF`                         | Enable ccache to speed up compilation |
| `DMRG_ENABLE_COVERAGE`            | `OFF`                         | Enable test coverage |
| `DMRG_BUILD_EXAMPLES`             | `OFF`                         | Build examples |
| `DMRG_PRINT_INFO`                 | `OFF`                         | Print information during cmake configure |
| `DMRG_PRINT_CHECKS`               | `OFF`                         | Print more information during cmake configure |
| `DMRG_DEPS_INSTALL_DIR`           | `CMAKE_INSTALL_PREFIX`        | Install directory for dependenciesx |
| `DMRG_DEPS_BUILD_DIR`             | `CMAKE_BINARY_DIR/dmrg-build` | Build directory for dependencies|
| `DMRG_PREFIX_ADD_PKGNAME`         | `OFF`                         | Install dependencies into CMAKE_INSTALL_PREFIX/<PackageName> |
| `BUILD_SHARED_LIBS`               | `OFF`                         | Link dependencies with static or shared libraries    |
| `CMAKE_INSTALL_PREFIX`            | None                          | Install directory for `h5pp` and dependencies  |
| `CONAN_PREFIX`                    | None                          | conan install directory  |

In addition, variables such
as [`<PackageName>_ROOT`](https://cmake.org/cmake/help/latest/variable/PackageName_ROOT.html)
and [`<PackageName>_DIR`](https://cmake.org/cmake/help/latest/command/find_package.html) can be set to help guide
CMake's `find_package` calls:


---

# Author

David Aceituno
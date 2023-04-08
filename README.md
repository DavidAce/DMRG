[![Ubuntu 22.04](https://github.com/DavidAce/DMRG/actions/workflows/ubuntu-22.04.yml/badge.svg)](https://github.com/DavidAce/DMRG/actions/workflows/ubuntu-22.04.yml)
[![codecov](https://codecov.io/gh/DavidAce/DMRG/branch/master/graph/badge.svg?token=9YE72CJ522)](https://codecov.io/gh/DavidAce/DMRG)
# DMRG++

---

# Introduction

[Density matrix renormalization group](https://en.wikipedia.org/wiki/Density_matrix_renormalization_group) (DMRG) is a
variational technique used to study quantum systems. DMRG works by iteratively optimizing a trial wave function in the
form of a [Matrix Product State](https://en.wikipedia.org/wiki/Matrix_product_states) (MPS), until obtaining an
eigenstate of the system with high precision.

DMRG++ includes 4 algorithms for finding eigenstates of 1D quantum spin chains:

- ***x*DMRG:** *Excited state* DMRG. Finds highly excited eigenstates of finite systems.
- ***f*DMRG:** *finite* DMRG. Finds the groundstate of finite systems.
- ***i*DMRG:** *infinite* DMRG. Finds the groundstate of infinite translationally invariant systems.
- ***i*TEBD:** *Imaginary Time Evolving Block Decimation*. Finds the ground state of infinite translationally
  invariant systems from a quench in imaginary time.

One additional algorithm is included to study the dynamics in the Many-body Localized phase:

- ***f*LBIT:** *Finite* l-BIT. Time evolution on a finite system in the basis of local integrals of motion (l-bits).


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
`DMRG++` takes a custom input file from command-line arguments, e.g. `./DMRG++ -c path/to/file.cfg`.  

## Output Data File

After execution the results are stored a binary file in HDF5 format. Its location is specified in the configuration
file `input/input.cfg`. By default, the output file is `output/output.h5`, which will contain values like the final
energies, entanglement entropies, and optionally the final state in MPS form.

To view the data you can use any hdf5-viewer, such as HDFCompass or HDFViewer.

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


## Requirements

The following software is required to build the project:

- C++17 compiler
- CMake version >= 3.24 to use conan as a CMake Dependency Provider. Otherwise, 3.19 is sufficient.

In addition, conan version 1.59 or higher is recommended for dependency installation. 
When using conan, you will need:

`conan remote add conan-dmrg https://neumann.theophys.kth.se.org/artifactory/api/conan/conan-dmrg`

to obtain arpack++ (see dependencies below).

## Quick start

- `git clone git@github.com:DavidAce/DMRG.git` and `cd DMRG`
- Configure `cmake --preset <preset>` (see available presets with `cmake --list-presets`)
- Build with `cmake --build --preset <preset>`
- Modify `input/input.config` to configure a simulation.
- Run with `./build/<preset>/DMRG++ -c input/input.cfg`.
- Find generated data in `output/output.h5`.

Some presets, with `conan` in their name, can use the 
[CMake Dependency Provider](https://cmake.org/cmake/help/latest/guide/using-dependencies/index.html#dependency-providers)
mechanism to let CMake call conan to install all the dependencies automatically.


## Dependencies

- Some BLAS, LAPACK and Lapacke implementation. Choose either [FlexiBLAS](https://www.mpi-magdeburg.mpg.de/projects/flexiblas) with reference Lapacke,  [Intel MKL](https://software.intel.com/en-us/mkl)
  or [OpenBLAS](https://github.com/xianyi/OpenBLAS). Use the [`BLA_VENDOR`](https://cmake.org/cmake/help/latest/module/FindBLAS.html) mechanism to guide CMake. to OpenBLAS can be built by conan. 
- [**Eigen**](http://eigen.tuxfamily.org) for tensor and matrix and linear algebra (tested with version >= 3.3).
- [**Arpack**](https://github.com/opencollab/arpack-ng) Eigenvalue solver based on Fortran. Note that this in turn
  requires LAPACK and BLAS libraries, both of which are included in OpenBLAS.
- [**Arpackpp**](https://github.com/m-reuter/eigsolver_properties) C++ frontend for Arpack.
- [**primme**](https://github.com/primme/primme) Eigenvalue solver. Complements Arpack.
- [**h5pp**](https://github.com/DavidAce/h5pp) a wrapper for HDF5. Includes [HDF5](https://www.hdfgroup.org/solutions/hdf5/), [spdlog](https://github.com/gabime/spdlog)
  and [fmt](https://github.com/fmtlib/fmt).
- [**ceres**](http://ceres-solver.org/) Optimization library with L-BFGS routines for unconstrained minimization.
- [**CLI11**](https://github.com/CLIUtils/CLI11) input argument parser
- [**Backward-cpp**](https://github.com/bombela/backward-cpp) pretty stack trace printer.


## Automatic Dependency Installation

The CMake flag `DMRG_PACKAGE_MANAGER` controls the automated behavior for finding or installing dependencies. It can
take one of these strings:

| Option               | Description                                                                                         |
|----------------------|-----------------------------------------------------------------------------------------------------|
| `find` **(default)** | Use CMake's `find_package`  to find dependencies. (Use this with the CMake Presets labeled `conan`) |
| `cmake¹`             | Use CMake to download and install dependencies during configure.                                    |
|                      |                                                                                                     |


## CMake Options

This project takes several flags in the form `cmake [-DOPTIONS=var] ../ `:

| Var                       | Default                       | Description                                                                         |
|---------------------------|-------------------------------|-------------------------------------------------------------------------------------|
| `DMRG_PACKAGE_MANAGER`    | `find`                        | Handle dependencies, `find` or `cmake`                                              |
| `DMRG_ENABLE_TBLIS`       | `OFF`                         | Use faster [tblis](https://github.com/devinamatthews/tblis) for tensor contractions |
| `DMRG_ENABLE_TESTS`       | `OFF`                         | Enable unit testing with ctest                                                      |
| `DMRG_ENABLE_DOCS`        | `OFF`                         | Build documentation                                                                 |
| `DMRG_ENABLE_COVERAGE`    | `OFF`                         | Enable test coverage                                                                |
| `DMRG_BUILD_EXAMPLES`     | `OFF`                         | Build examples                                                                      |
| `DMRG_BUILD_TOOLS`        | `OFF`                         | Build additional binaries under `./tools` for postprocessing (e.g. dmrg-meld)       |
| `DMRG_DEPS_INSTALL_DIR`   | `CMAKE_INSTALL_PREFIX`        | Install directory for dependencies                                                  |
| `DMRG_DEPS_BUILD_DIR`     | `CMAKE_BINARY_DIR/dmrg-build` | Build directory for dependencies                                                    |
| `DMRG_PREFIX_ADD_PKGNAME` | `OFF`                         | Install dependencies into CMAKE_INSTALL_PREFIX/<PackageName>                        |
| `DMRG_CMAKE_DEBUG`        | `OFF`                         | Extra information during CMake configuration                                        |
| `EIGEN_USE_THREADS`       | `ON`                          | Use STL threads to parallelize Eigen::Tensor (honors OMP_NUM_THREADS at runtime)    |
| `COMPILER_ENABLE_ASAN`    | `OFF`                         | Enable runtime address sanitizer -fsanitize=address                                 |
| `COMPILER_ENABLE_USAN`    | `OFF`                         | Enable undefined behavior sanitizer -fsanitize=undefined                            |
| `COMPILER_ENABLE_LTO`     | `OFF`                         | Enable link time optimization                                                       |
| `COMPILER_ENABLE_PCH`     | `OFF`                         | Enable precompiled headers to speed up compilation                                  |
| `COMPILER_ENABLE_CCACHE`  | `OFF`                         | Enable ccache to speed up compilation                                               |


In addition, variables such
as [`<PackageName>_ROOT`](https://cmake.org/cmake/help/latest/variable/PackageName_ROOT.html)
and [`<PackageName>_DIR`](https://cmake.org/cmake/help/latest/command/find_package.html) can be set to help guide
CMake's `find_package` calls:


---

# Contact Information
For questions about DMRG++ email David Aceituno aceituno `<at>` kth.se, or create a new [issue](https://github.com/DavidAce/DMRG/issues) or [discussion](https://github.com/DavidAce/DMRG/discussion).
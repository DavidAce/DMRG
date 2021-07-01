# Installation
## Quick start
- Git clone and build with `./build.sh` (see options with `./build.sh --help`).
- Modify `input/input.config` to configure a simulation.
- Run with `./build/Release/DMRG++ -c input/input.cfg` after setting your options in `input/input.config`.
- Find generated data in `output/output.h5`.

## Minimum Requirements
The following software is required to build the project:

- C++17 compiler. Tested with:
    - GNU GCC version >= 7
    - Clang version >= 6.0.
- * *To build dependencies*: Fortran compiler, tested with gfortran version >= 7
- CMake version >= 3.15.


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
²**                     | Use the [Conan package manager](https://conan.io/) to download and install dependencies automatically. Disregards libraries elsewhere on your system  |

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


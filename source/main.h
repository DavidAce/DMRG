#ifndef DMRG_MAIN_H
#define DMRG_MAIN_H

/*! \mainpage


  \brief This program finds the ground state of a 1D quantum Ising chain in a transverse field using the DMRG algorithm.

  Link to the [working notes](https://github.com/DavidAce/Notebooks/blob/master/DMRG%2B%2B/DMRG%2B%2B.pdf) on the theoretical aspects of this implementation.

  \tableofcontents

  \section intro Description of DMRG++

  [Density matrix renormalization group](https://en.wikipedia.org/wiki/Density_matrix_renormalization_group) (DMRG) is a variational numerical technique to study the low-energy physics of many-body quantum systems.

  This algorithm constructs and minimizes trial wave functions, in the shape of [Matrix Product States](https://en.wikipedia.org/wiki/Matrix_product_state) (MPS), iteratively in order to find the ground state of one-dimensional quantum systems with high precision.

  This implementation is inspired by the notation and steps in these articles:

  > [Phase Diagram of the Anisotropic Spin-2 XXZ Model: Infinite-System Density Matrix Renormalization Group Study](https://arxiv.org/abs/1212.6255)<br>
  > by Kjäll, Zaletel, Mong, Bardarson, and Pollmann.
  > Physical Review B 87 (23): 235106. <br>

  > [Efficient Numerical Simulations Using Matrix-Product States](http://quantumtensor.pks.mpg.de/wp-content/uploads/2016/06/notes_1.pdf)<br>
  > by Frank Pollmann. <br>

  > [The density-matrix renormalization group in the age of matrix product states](https://arxiv.org/abs/1008.3477)<br>
  > by Ulrich Schollwöck. <br>

  > [The iTEBD algorithm beyond unitary evolution](https://doi.org/10.1103/PhysRevB.78.155117)<br>
  > by Orus & Vidal. Arxiv, 070201 <br>

  > [Infinite size density matrix renormalization group, revisited](http://arxiv.org/abs/0804.2509)<br>
  > by McCulloch <br>

  > [Obtaining Highly Excited Eigenstates of Many-Body Localized Hamiltonians by the Density Matrix Renormalization Group Approach](https://doi.org/10.1103/PhysRevLett.116.247204)<br>
  > by Khemani, V., Pollmann, F., & Sondhi, S. L. (2016).
  > Physical Review Letters, 116(24), 1–5


---
\section installation Installation
\subsection quickstart Quick start
Git clone or copy & extract the project into a folder of your choosing.
**Make sure there are no spaces in your path!**.
The project can be built with a single command from a unix terminal.
Simply launch the script `.\build.sh` found in the root folder to trigger a CMake build.

The script takes optional arguments, run `.\build.sh -h` to learn more.

The CMake build will check for dependencies and download them automatically if needed (see *Optional Requirements* below).
If the dependencies are found, the project is built and an executable is generated.

To run the executable, launch `.\run.sh`, or  `.\run.sh -h` to learn more.


**Alternatively**, if you intend to develop or study the source code, some IDE's with CMake support can self-configure from the file CMakeLists.txt found in the project root folder. This
is perhaps an even simpler approach. Recommended: [CLion](https://www.jetbrains.com/clion/download) or [Visual Studio Code](https://code.visualstudio.com/) with the C++ and CMake Tools extensions.


\subsection minreqs Minimum Requirements
The following software is required to build the project:
 - C++ compiler with support for c++17 standard and libstdc++ standard library implementation  (version >= 7). Tested with two compilers:
    - GNU GCC versions 7 and 8 (these bundle libstdc++)
    - Clang version >= 5.0. (you need to manually install libstdc++ version >= 7, that comes bundled with gcc, for instance from `ppa:ubuntu-toolchain-r/test`)
 - CMake version >= 3.9. If you compile CMake from source, remember to enable `curl` (`./bootstrap --system-curl`).
 - *(For library compilation only) Fortran compiler, tested with gfortran-7 and gfortran-8.*

**Ubuntu** 17 or higher will have the versions required in the default repositories. For older distributions, use the ppa `ubuntu-toolchain-r/test` to get newer versions.

**Mac OSX** users are advised to use GNU GCC version 7 or 8 from homebrew. Install with `brew install gcc`. Clang from llvm 6.0 might work but you will have to link to GNU's `libstdc++.so` or `libstdc++.a` manually. The AppleClang compiler is not supported at all.


\subsection optionalreqs Optional Requirements
The compilation of DMRG++ requires several libraries. To meet the requirements, you have two options:

  1. **Automatic**: let cmake download and compile the libraries below from source into a local folder `libs/`. This is the default behavior if the library is not found on your system. Note that this does *NOT* make a system-wide install, so you can safely delete the `libs/` folder.
  2. **Manual**: install the libraries yourself with your favorite package manager (e.g. `apt` in ubuntu or `brew` in OSX). The build script will always attempt to find the libraries in your system first.

 The latter is recommended to avoid a lengthy compilation of these rather large libraries. If the compilation halts due to any library failing to compile or link, you can try installing/uninstalling that library from your package manager.

 #### List of libraries

 - **BLAS** and **LAPACK**. Required for Arpack. You can choose either [Intel MKL](https://software.intel.com/en-us/mkl) or [OpenBLAS](https://github.com/xianyi/OpenBLAS) (OpenBLAS requires Fortran to compile from source). If both MKL and OpenBLAS are found in the system, MKL is preferred.
 - [**Eigen**](http://eigen.tuxfamily.org) for tensor and matrix and linear algebra (tested with version >= 3.3).
 - [**Arpack**](https://github.com/opencollab/arpack-ng) Eigenvalue solver based on fortran. Note that this in turn requires LAPACK and BLAS libraries, both of which are included in OpenBLAS.
 - [**Arpackpp**](https://github.com/m-reuter/arpackpp) c++ frontend for Arpack.
 - [**Elemental**](http://libelemental.org/) for full diagonalization of matrices.
 - [**HDF5**](https://support.hdfgroup.org/HDF5/) for hdf5 binary output file support (tested with version >= 1.10).
 - [**GSL**](https://www.gnu.org/software/gsl/) for numerical integration (tested with version >= 2.4).

---



\section usage Usage

The executable `build/Release/DMRG++` can be run without input parameters. By default it will try to find `input/input.cfg` the file
where the simulation parameters are defined. You can modify these parameters, o create a new input file and pass its (full or relative) path as a command-line argument.

The results from the simulation are saved under `output/` as hdf5-files. To view the data you can use any hdf5-viewer, such as HDFCompass.

The script `Data_analysis/data_analysis.py` (in progress) shows how to analyze the simulation data stored in the hdf5 files. You need to install the python package
`h5py` from pip or conda to read files in the hdf5 format.

\section notation Notation

The *Vidal canonical form*, i.e. \f$\Gamma\Lambda\Gamma\Lambda\f$"..., is used throughout this code.
In code we denote

 - \f$\Gamma \rightarrow\f$ `G`.
 - \f$\Lambda \rightarrow\f$ `L`.

 \subsection convention Tensor index order convention.
 The tensor index order used here follows the convention:
 - physical indices first, from left to right or for MPO's, up to down.
 - other dimensions (like bond dimensions) next, from left to right.

\subsubsection example Example:
Consider for some position \f$n\f$ on the chain \f$\Gamma = \Gamma^{\sigma_n}_{a,b}\f$.
Here \f$\sigma_n \in [-1,1]\f$ is a particle with local (physical) dimension \f$d\f$ = 2, and \f$a,b\f$ are the remaining dimensions, in practice they are
bond dimensions of \f$\Lambda^{n-1}\f$ and \f$\Lambda^{n}\f$, respectively, which can be numbers \f$\in [1,\chi]\f$.

In diagrammatic tensor notation this is:
@verbatim
                 	    [d]          0
            G     =	[a]__|__[b] = 1__|__2
@endverbatim
where after the second equality the index order is shown. In code this corresponds to

\code{.cpp}
 using Scalar = std::complex<double>;
 long rank = 3
 Eigen::Tensor<Scalar,rank> G(d,a,b);
\endcode

An exception to this rule is theta: We have \f$\Theta^{\sigma_n,\sigma_{n+1}}_{\chi_a,\chi_b}\f$:

@verbatim
                 	           	[d] [d]                0   2
            Theta     =	[chia]___|___|___[chib] = 1 ___|___|___ 3
@endverbatim

which in code reads

\code{.cpp}
 using Scalar = std::complex<double>;
 long rank = 4
 Eigen::Tensor<Scalar,rank> theta(d,d,chia,chib);
\endcode

This index order doesn't follow the convention above because it simplifies the Schmidt-decomposition, where
the indices are merged in pairs ([0,1],[2,3]) to form a matrix on which to perform a singular value decomposition (SVD).



\section details Details
 \author    David Aceituno
 \date      02-2018
 \copyright MPL2

*/


#endif //DMRG_MAIN_H

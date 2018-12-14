#ifndef DMRG_MAIN_H
#define DMRG_MAIN_H

/*! \mainpage


  \brief This program finds the ground state of a 1D quantum system using the DMRG algorithm.

  [Working notes](https://github.com/DavidAce/Notebooks/blob/master/DMRG%2B%2B/DMRG%2B%2B.pdf) on the theoretical aspects of this implementation.

  \tableofcontents

  \section intro DMRG++

 [Density matrix renormalization group](https://en.wikipedia.org/wiki/Density_matrix_renormalization_group) (DMRG) is a variational numerical technique used to study many-body
  quantum systems. It works by optimizing a trial wave function in the form of a [Matrix Product States](https://en.wikipedia.org/wiki/Matrix_product_states) (MPS), to find either the
  groundstate or an eigenstate of a 1D quantum system with high precision. DMRG++ includes 4 different MPS-based algorithms for 1D systems:

  - **iDMRG:** *Infinite* DMRG. Finds the groundstate of infinite and translationally invariant systems.
  - **fDMRG:** *Finite* DMRG. Finds the groundstate of finite systems, not necessarily translationally invariant.
  - **xDMRG:** *Excited state* DMRG. Finds highly excited (mid-spectrum) eigenstates of finite systems.
  - **iTEBD:** *Imaginary Time Evolving Block Decimation*. Finds the ground state of infinite and translationally invariant systems using unitary operators that perform imaginary time evolution.

The program is controlled through an input file, whose path (full or relative to the binary) is given as input argument in the command line. See the [Installation](#installation) section below.

Included here are two 1D models of spin chains, the *Quantum Ising model with transverse field* and the *Self-dual quantum Ising model with random couplings and random fields*. The choice of model
is done in the input configuration file.

This implementation is inspired by the notation and steps in these articles:

  > [Phase Diagram of the Anisotropic Spin-2 XXZ Model: Infinite-System Density Matrix Renormalization Group Study](https://arxiv.org/abs/1212.6255)<br>
  > by Kjäll, Zaletel, Mong, Bardarson, and Pollmann, 2012 <br>

  > [Efficient Numerical Simulations Using Matrix-Product States](http://quantumtensor.pks.mpg.de/wp-content/uploads/2016/06/notes_1.pdf)<br>
  > by F. Pollmann. <br>

  > [The density-matrix renormalization group in the age of matrix product states](https://arxiv.org/abs/1008.3477)<br>
  > by U. Schollwöck, 2008 <br>

  > [The iTEBD algorithm beyond unitary evolution](https://doi.org/10.1103/PhysRevB.78.155117)<br>
  > by Orus & Vidal, 2008<br>

  > [Infinite size density matrix renormalization group, revisited](http://arxiv.org/abs/0804.2509)<br>
  > by McCulloch, 2008 <br>

  > [Obtaining Highly Excited Eigenstates of Many-Body Localized Hamiltonians by the Density Matrix Renormalization Group Approach](https://doi.org/10.1103/PhysRevLett.116.247204)<br>
  > by V. Khemani,  F. Pollmann  & S. L. Sondhi, 2016


\subsection notes Notes (in construction)
Go to the [working notes](https://github.com/DavidAce/Notebooks/blob/master/DMRG%2B%2B/DMRG%2B%2B.pdf) to learn more about the theoretical aspects of this implementation.


---
\section installation Installation
\subsection quickstart Quick start
- Git clone and build with `./build.sh`.
- Run with `./run.sh` after setting your options in `./input/input.cfg`.
- Profit from `./output/output.h5`.


\subsection requirements Requirements
The following software is required to build the project:
 - C++ compiler with support for c++17 standard and libstdc++ standard library implementation  (version >= 7). Tested with two compilers:
    - GNU GCC versions 7 and 8 (these bundle libstdc++)
    - Clang version >= 5.0. (you need to manually install libstdc++ version >= 7, that comes bundled with gcc, for instance from `ppa:ubuntu-toolchain-r/test`)
 - CMake version >= 3.9. If you compile CMake from source, remember to enable `curl` (`./bootstrap --system-curl`).
 - *For library compilation only:* Fortran compiler, tested with gfortran-7 and gfortran-8.

Ubuntu 17 or higher will have the versions required in the default repositories. For older distributions, use the ppa `ubuntu-toolchain-r/test` to get newer versions.

Mac OSX users are advised to use GNU GCC version 7 or 8 from homebrew. Install with `brew install gcc`. Clang from llvm 6.0 might work but you will have to link to GNU's `libstdc++.so` or `libstdc++.a` manually. The AppleClang compiler is not supported at all.


\subsubsection optional-requirements Optional Requirements
The compilation of DMRG++ requires several libraries. To meet the requirements, you have two options:

  1. **Automatic**: let cmake download and compile the libraries below from source into a local folder `libs/`. This is the default behavior if the library is not found on your system. Note that this does *NOT* make a system-wide install, so you can safely delete the `libs/` folder.
  2. **Manual**: install the libraries yourself with your favorite package manager (e.g. `apt` in ubuntu or `brew` in OSX). The build script will always attempt to find the libraries in your system first.

 The latter is recommended to avoid a lengthy compilation of these rather large libraries. If the compilation halts due to any library failing to compile or link, you can try installing/uninstalling that library from your package manager.

 #### List of libraries

 - **BLAS** and **LAPACK**. Required for Arpack. You can choose either [Intel MKL](https://software.intel.com/en-us/mkl) or [OpenBLAS](https://github.com/xianyi/OpenBLAS). If not found, OpenBLAS is downloaded automatically. Note that OpenBLAS requires Fortran to compile from source. If both MKL and OpenBLAS are found in the system, MKL is preferred.
 - [**Eigen**](http://eigen.tuxfamily.org) for tensor and matrix and linear algebra (tested with version >= 3.3).
 - [**Arpack**](https://github.com/opencollab/arpack-ng) Eigenvalue solver based on Fortran. Note that this in turn requires LAPACK and BLAS libraries, both of which are included in OpenBLAS.
 - [**Arpackpp**](https://github.com/m-reuter/eigsolver_properties) C++ frontend for Arpack.
 - [**HDF5**](https://support.hdfgroup.org/HDF5/) for hdf5 binary output file support (tested with version >= 1.10).
 - [**GSL**](https://www.gnu.org/software/gsl/) for numerical integration (tested with version >= 2.4).


\subsection build Build

Git clone or copy & extract the project into a folder of your choosing. Make sure there are no spaces in your path!

The project can be built with a single command from a unix terminal.
Simply launch the script `.\build.sh` found in the root folder to trigger a CMake build.

The script takes optional arguments, run `.\build.sh -h` to learn more.

CMake will check for dependencies in the host system. If not found, it will download and install these automatically to a folder `libs` in the project root (see *Optional Requirements* below).
If the dependencies are successfully found or installed, the project is built and an executable is generated in `./build/Release/DMRG++`.

Alternatively, if you intend to develop or study the source code, some IDE's with CMake support can self-configure from the file CMakeLists.txt found in the project root folder. This
is perhaps an even simpler approach. Recommended: [CLion](https://www.jetbrains.com/clion/download) or [Visual Studio Code](https://code.visualstudio.com/) with the C++ and CMake Tools extensions.

\subsection execution Execution
To run the executable, launch `.\run.sh`, or `.\run.sh -h` to learn more.
If no configuration file is given as argument to `run.sh`, the executable will look for a configuration file located in `./input/input.cfg`.

\subsection configuration-file Configuration file
The default configuration file in `./input/input.cfg` contains run-time instructions for DMRG++. Here you can choose the type of simulation, type of model, model parameters,
system size, precision as well as settings for profiling, data storage and console verbosity. Read the comments in the file to learn more.


\subsection output-file Output file
After execution the results are stored a binary file in HDF5 format. Its location is specified in the configuration file `./input/input.cfg`.
By default this should be in `.output/output.h5`. This file will contain values like the final energies, entanglement entropies, entanglement spectrums and
optionally the final state in MPS form.

To view the data you can use any hdf5-viewer, such as HDFCompass.

The script `analysis/data_analysis.py` (in progress) shows how to analyze the simulation data stored in the hdf5 files. You need to install the python package
`h5py` from pip or conda to read files in the HDF5 format.



---


\section notation Notation

The *Vidal canonical form*, i.e. ...\f$\Gamma\Lambda\Gamma\Lambda\f$..., is the underlying data structure for MPS throughout this code.
In code we denote

 - \f$\Gamma \rightarrow\f$ `G`
 - \f$\Lambda \rightarrow\f$ `L`

There are methods to obtain the MPS in *mixed canonical form* as well, i.e. \f$AAA...AACBB...BBB\f$,
where \f$A\f$'s are left unitary, \f$B\f$'s are right-unitary and \f$C\f$ is a (diagonal) bond-matrix.
The \f$A\f$'s and \f$B\f$'s are obtained from the Vidal canonical form by simple contraction:

 - \f$A = \Lambda \Gamma\f$
 - \f$B = \Gamma \Lambda\f$



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
 Eigen::Tensor<Scalar,rank> theta(d,chia,d,chib);
\endcode

This index order doesn't follow the convention above because it simplifies the Schmidt-decomposition, where
the indices are merged in pairs ([0,1],[2,3]) to form a matrix on which to perform a singular value decomposition (SVD).



\section details Details
 \author    David Aceituno
 \date      02-2018
 \copyright MPL2

*/


#endif //DMRG_MAIN_H

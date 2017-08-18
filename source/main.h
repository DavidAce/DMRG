#ifndef FINITE_DMRG_EIGEN_MAIN_H_H
#define FINITE_DMRG_EIGEN_MAIN_H_H

/*! \mainpage
 * \brief This program finds the ground state of a 1D quantum Ising chain in a transverse field using the DMRG algorithm.

  # DMRG (in development)
  [Density matrix renormalization group](https://en.wikipedia.org/wiki/Density_matrix_renormalization_group) (DMRG) is a variational numerical technique to study the low-energy physics of many-body quantum systems.

  This algorithm constructs and minimizes trial wave functions, in the shape of [Matrix Product States](https://en.wikipedia.org/wiki/Matrix_product_state) (MPS), iteratively in order to find the ground state of one-dimensional quantum systems with high precision.

  This implementation loosely follows the steps outlined in:

  > [Phase Diagram of the Anisotropic Spin-2 XXZ Model: Infinite-System Density Matrix Renormalization Group Study.](https://arxiv.org/abs/1212.6255)<br>
  > by Kjäll, Zaletel, Mong, Bardarson, and Pollmann. Physical Review B 87 (23): 235106. <br>

  > [Efficient Numerical Simulations Using Matrix-Product States](http://quantumtensor.pks.mpg.de/wp-content/uploads/2016/06/notes_1.pdf)<br>
  > by Frank Pollmann. <br>

  > [The density-matrix renormalization group in the age of matrix product states](https://arxiv.org/abs/1008.3477)<br>
  > by Ulrich Schollwöck. <br>



  ## Usage
  This code lacks an API or command-line parameters. As such, details of execution have to be
  set in source code, specifically one can edit model and simulation parameters in `source/n_model.h` and `source/n_settings.h`.

  The files `source/class_algorithms.cpp` and  `source/class_algorithms.h` contain routines for the infinite-DMRG,
  finite-DMRG and infinite-TEBD that can be called from `main.cpp`. The algorithms should ideally give similar
  ground-state energies, which is a good sanity check.

 ## Notation

 The *Vidal canonical form*, i.e. \f$\Gamma\Lambda\Gamma\Lambda\f$"..., is used throughout this code.
 In code we denote

 - \f$\Gamma \rightarrow\f$ `G`.
 - \f$\Lambda \rightarrow\f$ `L`.

 ## Tensor index order convention.
 The tensor index order used here follows the convention:
 - physical indices first, from left to right or for MPO's, up to down.
 - other dimensions (like bond dimensions) next, from left to right.

 #### Example:
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
 Textra::Tensor3 G(d,a,b);
\endcode

Similarly, we have for \f$\Theta^{\sigma_n,\sigma_{n+1}}_{\chi_a,\chi_b}\f$:

@verbatim
                 	           	[d] [d]                0   1
            Theta     =	[chia]___|___|___[chib] = 2 ___|___|___ 3
@endverbatim

which in code reads

\code{.cpp}
 Textra::Tensor4 G(d,d,chia,chib);
\endcode



 ## Requirements
 The following software is required:
 - [g++-7 for -std=c++17 support](https://launchpad.net/~ubuntu-toolchain-r/+archive/ubuntu/test) for newer features. To install in Ubuntu:
    \code
     sudo add-apt-repository ppa:ubuntu-toolchain-r/test
     sudo apt update
     sudo apt install gcc-7 g++-7
    \endcode
 - [Eigen](http://eigen.tuxfamily.org) for tensor support and SVD decomposition. **Included, no action needed**
 - [Spectra](https://spectralib.org/)  for eigenvalue-solver. **Included, no action needed**
 - [HDF5](https://support.hdfgroup.org/HDF5/) for hdf5-file storage support. To install in Ubuntu:
     \code
     sudo apt install libhdf5-dev libhdf5-10 libhdf5-cpp-11
     \endcode
 - [GSL](https://www.gnu.org/software/gsl/) for numerical integration. To install in Ubuntu:
    \code
    sudo apt install libgsl-dev libgsl2
    \endcode


# Details
 \author    David Aceituno
 \date      07-2017
 \copyright MPL2

*/


#endif //FINITE_DMRG_EIGEN_MAIN_H_H

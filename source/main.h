#ifndef DMRG_MAIN_H
#define DMRG_MAIN_H

/*! \mainpage


  \brief This program finds the ground state of a 1D quantum system using the DMRG algorithm.

  [Working notes](https://github.com/DavidAce/Notebooks/blob/master/DMRG%2B%2B/DMRG%2B%2B.pdf) on the theoretical aspects of this implementation.

  \tableofcontents

  \section intro DMRG++ (documentation in progress)

[Density matrix renormalization group](https://en.wikipedia.org/wiki/Density_matrix_renormalization_group) (DMRG) is a variational numerical technique used to
study many-body quantum systems. It works by optimizing a trial wave function in the form of a [Matrix Product
States](https://en.wikipedia.org/wiki/Matrix_product_states) (MPS), to find either the groundstate or an eigenstate of a 1D quantum system with high precision.
DMRG++ includes 4 different MPS-based algorithms for 1D systems:

  - **iDMRG:** *infinite* DMRG. Finds the groundstate of infinite and translationally invariant systems.
  - **fDMRG:** *finite* DMRG. Finds the groundstate of finite systems, not necessarily translationally invariant.
  - **xDMRG:** *Excited state* DMRG. Finds highly excited (mid-spectrum) eigenstates of finite systems.
  - **iTEBD:** *Imaginary Time Evolving Block Decimation*. Finds the ground state of infinite and translationally invariant systems using unitary operators that
perform imaginary time evolution.

The program is controlled through an input file, whose path (full or relative to the binary) is given as input argument in the command line. See the
[Installation](#installation) section below.

Included here are two 1D models of spin chains, the *Quantum Ising model with transverse field* and the *Self-dual quantum Ising model with random couplings and
random fields*. The choice of model is done in the input configuration file.


\subsection notes Notes (in construction)
Go to the [working notes](https://github.com/DavidAce/Notebooks/blob/master/DMRG%2B%2B/DMRG%2B%2B.pdf) to learn more about the theoretical aspects of this
implementation.


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

\subsection reference Reference
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

  > [Obtaining Highly Excited Eigenstates of Many-Body Localized Hamiltonians by the Density Matrix Renormalization Group
Approach](https://doi.org/10.1103/PhysRevLett.116.247204)<br> > by V. Khemani,  F. Pollmann  & S. L. Sondhi, 2016



\section details Details
 \author    David Aceituno
 \date      06-2019
 \copyright MPL2

*/

#endif // DMRG_MAIN_H

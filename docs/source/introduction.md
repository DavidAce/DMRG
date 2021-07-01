# Introduction

[Density matrix renormalization group](https://en.wikipedia.org/wiki/Density_matrix_renormalization_group) (DMRG) is a variational technique used to study 1D quantum systems. It works by iteratively optimizing a trial wave function in the form of a [Matrix Product State](https://en.wikipedia.org/wiki/Matrix_product_states) (MPS), until obtaining an eigenstate of the system with high precision. DMRG++ includes 4 different algorithms:

- ***i*DMRG:** *infinite* DMRG. Finds the groundstate of infinite and translationally invariant systems.
- ***f*DMRG:** *finite* DMRG. Finds the groundstate of finite systems, not necessarily translationally invariant.
- ***x*DMRG:** *Excited state* DMRG. Finds highly excited (mid-spectrum) eigenstates of finite systems.
- ***f*LBIT:** *Finite* l-BIT. Time evolution on a finite system in the basis of local integrals of motion (the
  l-bits) of an MBL phase.
- ***i*TEBD:** *Imaginary Time Evolving Block Decimation*. Finds the ground state of infinite and translationally
  invariant systems using unitary operators that perform imaginary time evolution.

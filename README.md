 
 ### Documentation
 See the [Documentation](https://davidace.github.io/DMRG/) page generated by [Doxygen](www.doxygen.org).


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


---
## Quick Start

#### From command line
To build project, enter:
```
./build.sh
```

To run executable, enter:

```
./run.sh
```

#### From IDE
Some IDE's with CMake support can self-configure from the file CMakeLists.txt found in the project root folder. This
is by far the easiest approach. Recommended: [CLion](https://www.jetbrains.com/clion/download) or [Visual Studio Code](https://code.visualstudio.com/) with C++ and CMake Tools extensions.



#### Requirements
 Please install the following software before building the project.
 * g++ compiler with -std=c++17 support  (tested with g++ version >= 6)
 * CMake (tested with version >= 3.7)
 * [~~Spectra~~](https://spectralib.org/) ~~for eigenvalue-solver.~~ **Included, no action needed**
 
 
 The package manager [Hunter](https://github.com/ruslo/hunter) is included to ease the building process.
 During the first build, the dependencies listed in CMakeLists.txt will be downloaded and installed by
 [Hunter](https://github.com/ruslo/hunter) automatically on any platform (Linux/OSX/Win).
 
 The following software is installed by [Hunter](https://github.com/ruslo/hunter):   
 * [Eigen](http://eigen.tuxfamily.org) for tensor support and SVD decomposition.
 * [HDF5](https://support.hdfgroup.org/HDF5/) for hdf5-file storage support. 
 * [GSL](https://www.gnu.org/software/gsl/) for numerical integration.
 
 The default installation folder in Linux is `~/.hunter`.

---
 
## Usage
This code lacks an API or command-line parameters. As such, details of execution have to be
set in source code, specifically one can edit model and simulation parameters in `source/n_model.h` and `source/n_settings.h`.

The files `source/class_algorithms.cpp` and  `source/class_algorithms.h` contain routines for the infinite-DMRG,
finite-DMRG and infinite-TEBD that can be called from `main.cpp`. The algorithms should ideally give similar 
ground state energies, which is a good sanity check.

 ## Tensor index order convention.
 The tensor index order used here follows the convention:
 - physical indices first, from left to right or for MPO's, up to down.
 - other dimensions (like bond dimensions) next, from left to right.

 #### Example:
Consider a rank-3 tensor `G` with dimensions `d`, `a` and `b`. In diagrammatic tensor notation this is:
```
                 	    [d]          0
            G     =	[a]__|__[b] = 1__|__2
```

where after the second equality the index order is shown. In code this corresponds to

```
 Textra::Tensor3d G(d,a,b);
```

Similarly, we have for a rank-4 tensor `Theta`:

```
                                    [d] [d]                0   1
                Theta   =   [chiL]___|___|___[chiR] = 2 ___|___|___ 3
```

which in code reads

```
 Textra::Tensor4d G(d,d,chiL,chiR);
```




## License
Open source under [MPL2](https://www.mozilla.org/MPL/2.0/).

## Author
David Aceituno
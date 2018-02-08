 [![Build Status](https://travis-ci.org/DavidAce/DMRG.svg?branch=master)](https://travis-ci.org/DavidAce/DMRG)
 
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
## Installation
### Quick start
Git clone or copy & extract the project into a folder of your choosing.
**Make sure there are no spaces in the output_folder!**.
The project can be built with a single command from a unix terminal. 
Simply launch the script `.\build.sh` found in the root folder, containing
```
        #!/bin/sh
        cmake -E make_directory build/Release
        cd build/Release
        cmake -Bbuild/Release --build build -config Release ../../
        make
```
This script will create subdirectories and use CMake to check for dependencies and download them automatically if needed (see *Optional Requirements* below).
If the dependencies are found, the project is built and an executable is generated.

To run executable, launch `.\run.sh`, containing

```
#!/bin/sh 
./build/Release/DMRG++
```


**Alternatively** some IDE's with CMake support can self-configure from the file CMakeLists.txt found in the project root folder. This
is perhaps an even simpler approach. Recommended: [CLion](https://www.jetbrains.com/clion/download) or [Visual Studio Code](https://code.visualstudio.com/) with C++ and CMake Tools extensions.


### Minimum Requirements
The following software is required to build the project:
 - C++ compiler with support for c++17 standard and libstdc++ standard library implementation. Tested with both
    - GNU GCC version >= 7
    - Clang version >= 5.0
 - Fortran compiler. Tested with GNU GFORTRAN version >= 4. This is only needed to install libraries from source, see below.
 - CMake version >= 3.8. If you compile CMake from source, remember to enable `curl` (`./bootstrap --system-curl`). 
 
The AppleClang compiler is not supported.

### Optional Requirements
 You can chose to **either** 
  - let the build script compile the libraries below from source into `libs/`. This will happen automatically if the library is not found on your system. Note that this does *NOT* make a system-wide install, so you can safely delete the `libs/` folder.
  - (*recommended*) install the libraries yourself with your favorite package manager (e.g. `apt` in ubuntu or `brew` in OSX). The build script will always attempt to find the libraries in your system first.
 
 The latter is recommended to avoid a lengthy compilation of these rather large libraries. 
 
 #### Libraries
 
 - [Eigen](http://eigen.tuxfamily.org) for tensor support and SVD decomposition (tested with version >= 3.3).
 - [Arpack](https://github.com/opencollab/arpack-ng) Eigenvalue solver based on fortran. Note that this in turn requires LAPACK and BLAS API.
 - [Arpackpp](https://github.com/m-reuter/arpackpp) C++ frontend for ARPACK.
 - [Spectra](https://spectralib.org/) Another eigenvalue solver. Header only.
 - [HDF5](https://support.hdfgroup.org/HDF5/) for hdf5-file env_storage support (tested with version >= 1.8).
 - [GSL](https://www.gnu.org/software/gsl/) for numerical integration (tested with version >= 2.4).



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

For more information about tensors see the documentation for the [Eigen Tensor Module](https://bitbucket.org/eigen/eigen/src/e8005fc30c6956e3f413a8d7aa2dd6395f330ffe/unsupported/Eigen/CXX11/src/Tensor/README.md?at=default&fileviewer=file-view-default).


## License
Open source under [MPL2](https://www.mozilla.org/MPL/2.0/).

## Author
David Aceituno
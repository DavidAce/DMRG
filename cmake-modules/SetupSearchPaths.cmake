# Let cmake find our Find<package>.cmake modules
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake-modules)
if (CMAKE_SIZEOF_VOID_P EQUAL 8 OR CMAKE_GENERATOR MATCHES "64")
    set(FIND_LIBRARY_USE_LIB64_PATHS ON)
elseif (CMAKE_SIZEOF_VOID_P EQUAL 4)
    set(FIND_LIBRARY_USE_LIB32_PATHS ON)
endif ()

set(CMAKE_FIND_DEBUG_MODE OFF)

# Manually installed libraries are installed to CMAKE_INSTALL_PREFIX
list(APPEND CMAKE_PREFIX_PATH ${CMAKE_INSTALL_PREFIX})

# Paths to search for conda libraries
# These paths should only be searched when H5PP_PREFER_CONDA_LIBS = ON
if (DMRG_PREFER_CONDA_LIBS)
    list(APPEND DMRG_CONDA_CANDIDATE_PATHS
            ${CONDA_PREFIX}
            $ENV{CONDA_PREFIX}
            $ENV{HOME}/anaconda3/envs/dmrg
            $ENV{HOME}/anaconda/envs/dmrg
            $ENV{HOME}/miniconda3/envs/dmrg
            $ENV{HOME}/miniconda/envs/dmrg
            $ENV{HOME}/.conda/envs/dmrg
            $ENV{HOME}/anaconda3
            $ENV{HOME}/anaconda
            $ENV{HOME}/miniconda3
            $ENV{HOME}/miniconda
            $ENV{HOME}/.conda
            )
endif ()

if (DMRG_DOWNLOAD_METHOD MATCHES "conan")
    # Paths to search for conan installation.
    list(APPEND DMRG_CONAN_CANDIDATE_PATHS
            ${CONAN_PREFIX}
            $ENV{CONAN_PREFIX}
            ${CONDA_PREFIX}
            $ENV{CONDA_PREFIX}
            $ENV{HOME}/anaconda3/envs/dmrg
            $ENV{HOME}/anaconda/envs/dmrg
            $ENV{HOME}/miniconda3/envs/dmrg
            $ENV{HOME}/miniconda/envs/dmrg
            $ENV{HOME}/.conda/envs/dmrg
            $ENV{HOME}/anaconda3
            $ENV{HOME}/anaconda
            $ENV{HOME}/miniconda3
            $ENV{HOME}/miniconda
            $ENV{HOME}/.conda
            )
endif ()

if (DMRG_DOWNLOAD_METHOD MATCHES "find")
    # If modules are loaded check those paths as well
    list(APPEND CMAKE_PREFIX_PATH
            $ENV{EBROOTIMKL}
            $ENV{EBROOTOPENBLAS}
            $ENV{EBROOTBLAS}
            $ENV{EBROOTLAPACK}
            $ENV{EBROOTARPACKMINNG}
            $ENV{EBROOTARPACKPLUSPLUS}
            $ENV{EBROOTGCC}
            $ENV{EBROOTOPENMPI}
            $ENV{EBROOTZLIB}
            $ENV{EBROOTGLOG}
            $ENV{EBROOTGFLAGS}
            $ENV{EBROOTCERES}
            $ENV{EBROOTSUITESPARSE}
            $ENV{EBROOTCXSPARSE}
            $ENV{EBROOTMETIS}
            $ENV{EBROOTCHOLMOD}
            $ENV{EBROOTCOLAMD}
            $ENV{EBROOTCCOLAMD}
            $ENV{EBROOTAMD}
            $ENV{EBROOTCAMD}
            $ENV{BLAS_DIR}
            $ENV{BLAS_ROOT}
            $ENV{ARPACKPP_DIR}
            )
endif ()
if(DMRG_DOWNLOAD_METHOD MATCHES "find")
    # Let cmake find our Find<package>.cmake modules
    list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

    if(DMRG_PREFER_CONDA_LIBS)
        list(APPEND CMAKE_PREFIX_PATH
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
        endif()
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

endif()
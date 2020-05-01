if(DMRG_PRINT_INFO)
    # Print CMake options
    message(STATUS  "|----------------\n"
            "-- | BUILD_SHARED_LIBS       : ${BUILD_SHARED_LIBS}\n"
            "-- | CMAKE_INSTALL_PREFIX    : ${CMAKE_INSTALL_PREFIX}\n"
            "-- | DMRG_BUILD_EXAMPLES     : ${DMRG_BUILD_EXAMPLES}\n"
            "-- | DMRG_ENABLE_TESTS       : ${DMRG_ENABLE_TESTS}\n"
            "-- | DMRG_ENABLE_OPENMP      : ${DMRG_ENABLE_OPENMP}\n"
            "-- | DMRG_ENABLE_MKL         : ${DMRG_ENABLE_MKL}\n"
            "-- | DMRG_ENABLE_LTO         : ${DMRG_ENABLE_LTO}\n"
            "-- | DMRG_DOWNLOAD_METHOD    : ${DMRG_DOWNLOAD_METHOD}\n"
            "-- | DMRG_PREFER_CONDA_LIBS  : ${DMRG_PREFER_CONDA_LIBS}\n"
            "-- | DMRG_PRINT_INFO         : ${DMRG_PRINT_INFO}\n")
endif ()

if(DMRG_PRINT_INFO)
    # Print CMake options
    message(STATUS  "|----------------\n"
            "-- | CMAKE_INSTALL_PREFIX   : ${CMAKE_INSTALL_PREFIX}\n"
            "-- | CMAKE_VERBOSE_MAKEFILE : ${CMAKE_VERBOSE_MAKEFILE}\n"
            "-- | BUILD_SHARED_LIBS      : ${BUILD_SHARED_LIBS}\n"
            "-- | DMRG_ENABLE_THREADS    : ${DMRG_ENABLE_THREADS}\n"
            "-- | DMRG_ENABLE_OPENMP     : ${DMRG_ENABLE_OPENMP}\n"
            "-- | DMRG_ENABLE_MKL        : ${DMRG_ENABLE_MKL}\n"
            "-- | DMRG_ENABLE_TESTS      : ${DMRG_ENABLE_TESTS}\n"
            "-- | DMRG_ENABLE_ASAN       : ${DMRG_ENABLE_ASAN}\n"
            "-- | DMRG_ENABLE_USAN       : ${DMRG_ENABLE_USAN}\n"
            "-- | DMRG_ENABLE_LTO        : ${DMRG_ENABLE_LTO}\n"
            "-- | DMRG_ENABLE_PCH        : ${DMRG_ENABLE_PCH}\n"
            "-- | DMRG_ENABLE_CCACHE     : ${DMRG_ENABLE_CCACHE}\n"
            "-- | DMRG_BUILD_EXAMPLES    : ${DMRG_BUILD_EXAMPLES}\n"
            "-- | DMRG_DOWNLOAD_METHOD   : ${DMRG_DOWNLOAD_METHOD}\n"
            "-- | DMRG_PRINT_INFO        : ${DMRG_PRINT_INFO}\n"
            "-- | DMRG_PRINT_CHECKS      : ${DMRG_PRINT_CHECKS}\n"
            "-- | DMRG_DEPS_IN_SUBDIR    : ${DMRG_DEPS_IN_SUBDIR}\n"
            "-- | DMRG_PREFER_CONDA_LIBS : ${DMRG_PREFER_CONDA_LIBS}\n")
endif ()

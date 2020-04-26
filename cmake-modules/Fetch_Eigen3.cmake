
if(NOT TARGET Eigen3::Eigen AND DMRG_DOWNLOAD_METHOD MATCHES "find|fetch|native")
    # We want to find our own Eigen3 to make sure we patch it properly
    find_package(Eigen3
        HINTS ${CMAKE_INSTALL_PREFIX}/Eigen3
        PATH_SUFFIXES Eigen3 eigen3 include/Eigen3 include/eigen3
        NO_DEFAULT_PATH)
    if(TARGET Eigen3::Eigen)
        message(STATUS "Found Eigen3: ${EIGEN3_INCLUDE_DIR}")
        target_include_directories(Eigen3::Eigen SYSTEM INTERFACE ${EIGEN3_INCLUDE_DIR})
    endif()
endif()

if(NOT TARGET Eigen3::Eigen AND DMRG_DOWNLOAD_METHOD MATCHES "none")
    message(FATAL_ERROR "Eigen3 has to be downloaded because we need a patched version, "
            "which isn't available through any systems package manager")
endif()

if(NOT TARGET Eigen3::Eigen AND DMRG_DOWNLOAD_METHOD MATCHES "fetch|native")
    message(STATUS "Eigen3 will be installed into ${CMAKE_INSTALL_PREFIX}")
    include(${PROJECT_SOURCE_DIR}/cmake-modules/BuildDependency.cmake)
    build_dependency(Eigen3 "${CMAKE_INSTALL_PREFIX}" "")
    find_package(Eigen3 3.3.7
            HINTS ${CMAKE_INSTALL_PREFIX}/Eigen3
            PATH_SUFFIXES Eigen3 eigen3 include/Eigen3 include/eigen3
            NO_CMAKE_PACKAGE_REGISTRY NO_DEFAULT_PATH)
    if(TARGET Eigen3::Eigen)
        message(STATUS "Eigen3 installed successfully")
        target_include_directories(Eigen3::Eigen SYSTEM INTERFACE ${EIGEN3_INCLUDE_DIR})
    else()
        message(WARNING "Eigen3 could not be installed")
    endif()
endif()


if(TARGET Eigen3::Eigen AND TARGET openmp::openmp)
    target_compile_definitions    (Eigen3::Eigen INTERFACE -DEIGEN_USE_THREADS)
endif()

if(TARGET Eigen3::Eigen AND TARGET blas::blas )
    set(EIGEN3_USING_BLAS ON)
    if(TARGET mkl::mkl)
        message(STATUS "Eigen3 will use MKL")
        target_compile_definitions    (Eigen3::Eigen INTERFACE -DEIGEN_USE_MKL_ALL)
        target_compile_definitions    (Eigen3::Eigen INTERFACE -DEIGEN_USE_LAPACKE_STRICT)
        target_link_libraries         (Eigen3::Eigen INTERFACE mkl::mkl)
    else ()
        message(STATUS "Eigen3 will use OpenBLAS")
        target_compile_definitions    (Eigen3::Eigen INTERFACE -DEIGEN_USE_BLAS)
        target_compile_definitions    (Eigen3::Eigen INTERFACE -DEIGEN_USE_LAPACKE_STRICT)
        target_link_libraries         (Eigen3::Eigen INTERFACE blas::blas)
    endif()


    # Use this flag if Ceres is giving you trouble!
    # For some reason it starts mixing aligned and hand-made aligned malloc and freeing them willy nilly
    # This flag forces its hand and avoids a segfault in some cases.
    #    target_compile_definitions(Eigen3::Eigen INTERFACE -DEIGEN_MALLOC_ALREADY_ALIGNED=0) # Finally something works!!!


endif()


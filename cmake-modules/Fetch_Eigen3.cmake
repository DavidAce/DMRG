
if(NOT TARGET Eigen3::Eigen AND DMRG_DOWNLOAD_METHOD STREQUAL "find")
    message(STATUS "DMRG needs a patched version of Eigen 3.3.7 and will downloaded despite DMRG_DOWNLOAD_METHOD=find")
endif()

if(NOT TARGET Eigen3::Eigen AND DMRG_DOWNLOAD_METHOD MATCHES "find|fetch")
    # We want to find our own Eigen3 to make sure we patch it properly
    find_package(Eigen3
        HINTS ${CMAKE_INSTALL_PREFIX}
        NO_DEFAULT_PATH) # IMPORTANT TO ONLY LOOK IN DMRG'S OWN PLACE
    if(TARGET Eigen3::Eigen)
        message(STATUS "Found Eigen3: ${EIGEN3_INCLUDE_DIR}")
        target_include_directories(Eigen3::Eigen SYSTEM INTERFACE ${EIGEN3_INCLUDE_DIR})
    endif()
endif()

if(NOT TARGET Eigen3::Eigen AND DMRG_DOWNLOAD_METHOD MATCHES "find|fetch")
    message(STATUS "Eigen3 will be installed into ${CMAKE_INSTALL_PREFIX}")
    include(${PROJECT_SOURCE_DIR}/cmake-modules/BuildDependency.cmake)
    build_dependency(Eigen3 "${CMAKE_INSTALL_PREFIX}" "")
    find_package(Eigen3 3.3.7
            HINTS ${CMAKE_INSTALL_PREFIX}
            NO_DEFAULT_PATH)
    if(TARGET Eigen3::Eigen)
        message(STATUS "Eigen3 installed successfully")
        target_include_directories(Eigen3::Eigen SYSTEM INTERFACE ${EIGEN3_INCLUDE_DIR})
    else()
        message(FATAL_ERROR "Eigen3 could not be installed")
    endif()
endif()


if(TARGET Eigen3::Eigen)
    message(STATUS "Applying special Eigen compile definitions")

    # AVX aligns 32 bytes (AVX512 aligns 64 bytes).
    # When running on Tetralith, with march=native, there can be alignment mismatch
    # in ceres which results in a segfault on free memory.
    # Something like "double free or corruption ..."
    #   * EIGEN_MAX_ALIGN_BYTES=16 works on Tetralith
    cmake_host_system_information(RESULT _host_name   QUERY HOSTNAME)
    if(${_host_name} MATCHES "etralith|riolith")
        #target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_MALLOC_ALREADY_ALIGNED=0) # May work to fix CERES segfaults!!!
        target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_MAX_ALIGN_BYTES=16)
    endif()

endif()


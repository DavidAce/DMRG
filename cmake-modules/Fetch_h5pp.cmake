if(NOT TARGET h5pp::h5pp AND DMRG_DOWNLOAD_METHOD MATCHES "find|fetch|native")
    set(EIGEN3_NO_CMAKE_PACKAGE_REGISTRY TRUE)
    set(SPDLOG_NO_CMAKE_PACKAGE_REGISTRY TRUE)
    find_package(h5pp 1.7.0
            HINTS ${CMAKE_INSTALL_PREFIX} ${CONDA_HINTS}
            PATHS ${CMAKE_INSTALL_PREFIX}
            PATH_SUFFIXES h5pp
            NO_CMAKE_PACKAGE_REGISTRY)
    if(h5pp_FOUND AND TARGET h5pp::h5pp)
        message(STATUS "Found h5pp")
    endif()
endif()

if(NOT TARGET h5pp::h5pp AND DMRG_DOWNLOAD_METHOD MATCHES "fetch|native")
    message(STATUS "h5pp will be installed into ${CMAKE_INSTALL_PREFIX}")
    include(${PROJECT_SOURCE_DIR}/cmake-modules/BuildDependency.cmake)
    if(TARGET Eigen3::Eigen)
        get_target_property(EIGEN3_INCLUDE_DIR Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
        list(APPEND H5PP_CMAKE_OPTIONS  -DEigen3_DIR:PATH=${CMAKE_INSTALL_PREFIX}/Eigen3)
        list(APPEND H5PP_CMAKE_OPTIONS  -DEIGEN3_INCLUDE_DIR:PATH=${CMAKE_INSTALL_PREFIX}/Eigen3/include/eigen3)
    endif()
    list(APPEND H5PP_CMAKE_OPTIONS  -DH5PP_DIRECTORY_HINTS:PATH=${CMAKE_INSTALL_PREFIX})
    list(APPEND H5PP_CMAKE_OPTIONS  -DH5PP_DOWNLOAD_METHOD:BOOL=${DMRG_DOWNLOAD_METHOD})
    list(APPEND H5PP_CMAKE_OPTIONS  -DH5PP_PREFER_CONDA_LIBS:BOOL=${DMRG_PREFER_CONDA_LIBS})
    build_dependency(h5pp "${CMAKE_INSTALL_PREFIX}" "${H5PP_CMAKE_OPTIONS}")

    find_package(h5pp 1.7.0 HINTS ${CMAKE_INSTALL_PREFIX} PATHS PATH_SUFFIXES h5pp NO_DEFAULT_PATH)
    if(h5pp_FOUND AND TARGET h5pp::h5pp)
        message(STATUS "h5pp installed successfully")
    else()
        message(FATAL_ERROR "h5pp could not be installed")
    endif()
endif()

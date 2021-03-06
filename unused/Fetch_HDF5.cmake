include(../cmake-modules/FindPackageHDF5.cmake)

if(NOT TARGET hdf5::hdf5)
    find_package_hdf5()
    if(TARGET hdf5::hdf5)
#        message(STATUS "hdf5 found")

    elseif (DOWNLOAD_MISSING)
        message(STATUS "HDF5 will be installed into ${CMAKE_INSTALL_PREFIX}")
        include(../cmake-modules/BuildDependency.cmake)
        build_dependency(hdf5 "")
        set(HDF5_ROOT ${CMAKE_BINARY_DIR}/h5pp-deps-install)
        find_package_hdf5()
        if(TARGET hdf5::hdf5)
            message(STATUS "hdf5 installed successfully: ${HDF5_BUILD_DIR} ${HDF5_CXX_INCLUDE_DIRS} ${HDF5_hdf5_LIBRARY}")
        else()
            message(STATUS "config_result: ${config_result}")
            message(STATUS "build_result: ${build_result}")
            message(FATAL_ERROR "hdf5 could not be downloaded.")
        endif()
    else()
        message("WARNING: Dependency HDF5 not found and DOWNLOAD_MISSING is OFF. Build may fail.")
    endif()

endif()


find_package(h5pp PATHS ${CMAKE_INSTALL_PREFIX}/h5pp $ENV{H5PP_DIR})

if(h5pp_FOUND AND TARGET h5pp::h5pp)
    message(STATUS "h5pp found")
elseif(DOWNLOAD_MISSING)
    message(STATUS "h5pp will be installed into ${CMAKE_INSTALL_PREFIX}")
    include(${PROJECT_SOURCE_DIR}/cmake-modules/BuildDependency.cmake)
    get_target_property(EIGEN3_INCLUDE_DIR Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
    list(APPEND H5PP_CMAKE_OPTIONS  -DEigen3_DIR:PATH=${CMAKE_INSTALL_PREFIX}/Eigen3)
    list(APPEND H5PP_CMAKE_OPTIONS  -DEIGEN3_INCLUDE_DIR:PATH=${CMAKE_INSTALL_PREFIX}/Eigen3/include/eigen3)
    build_dependency(h5pp "${H5PP_CMAKE_OPTIONS}")
    find_package(h5pp PATHS ${CMAKE_INSTALL_PREFIX}/h5pp $ENV{H5PP_DIR})
    if(h5pp_FOUND AND TARGET h5pp::h5pp)
        message(STATUS "h5pp installed successfully")
    else()
        message(STATUS "config_result: ${config_result}")
        message(STATUS "build_result: ${build_result}")
        message(FATAL_ERROR "h5pp could not be downloaded.")
    endif()
else()
    message(FATAL_ERROR "Dependency h5pp not found and DOWNLOAD_MISSING is OFF")
endif()

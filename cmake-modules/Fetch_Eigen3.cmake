# We want to find our own Eigen3 to make sure we patch it properly
find_package(Eigen3
        HINTS ${CMAKE_INSTALL_PREFIX}/Eigen3
        PATHS ${CMAKE_INSTALL_PREFIX}/Eigen3
        PATH_SUFFIXES Eigen3 eigen3 include/Eigen3 include/eigen3  NO_CMAKE_PACKAGE_REGISTRY NO_DEFAULT_PATH)

if(TARGET Eigen3::Eigen)
    message(STATUS "Eigen3 found")
    target_link_libraries(Eigen3::Eigen INTERFACE -lrt)
#    include(cmake-modules/PrintTargetProperties.cmake)
#    print_target_properties(Eigen3::Eigen)

elseif (DOWNLOAD_MISSING)
    message(STATUS "Eigen3 will be installed into ${CMAKE_INSTALL_PREFIX}")
    include(${PROJECT_SOURCE_DIR}/cmake-modules/BuildDependency.cmake)
    build_dependency(Eigen3 "")
    find_package(Eigen3 3.3.7
            HINTS ${CMAKE_INSTALL_PREFIX}/Eigen3
            PATHS ${CMAKE_INSTALL_PREFIX}/Eigen3
            PATH_SUFFIXES Eigen3 eigen3 include/Eigen3 include/eigen3  NO_CMAKE_PACKAGE_REGISTRY NO_DEFAULT_PATH)
    if(TARGET Eigen3::Eigen)
        message(STATUS "Eigen3 installed successfully")
        target_link_libraries(Eigen3::Eigen INTERFACE -lrt)
    else()
        message(STATUS "cfg_result: ${cfg_result}")
        message(STATUS "bld_result: ${bld_result}")
        message(FATAL_ERROR "Eigen3 could not be downloaded.")
    endif()

else()
    message(STATUS "Dependency Eigen3 not found and DOWNLOAD_MISSING is OFF")
endif()


if(TARGET Eigen3::Eigen AND TARGET blas )
    set(EIGEN3_USING_BLAS ON)
    if(TARGET mkl)
        message(STATUS "Eigen3 will use MKL")
        target_compile_definitions    (Eigen3::Eigen INTERFACE -DEIGEN_USE_MKL_ALL)
        target_compile_definitions    (Eigen3::Eigen INTERFACE -DEIGEN_USE_LAPACKE_STRICT)
        target_compile_definitions    (Eigen3::Eigen INTERFACE -DEIGEN_USE_MKL_ALL)
        target_include_directories(Eigen3::Eigen INTERFACE ${MKL_INCLUDE_DIR})
    else ()
        message(STATUS "Eigen3 will use BLAS and LAPACKE")
        target_compile_definitions    (Eigen3::Eigen INTERFACE -DEIGEN_USE_BLAS)
        target_compile_definitions    (Eigen3::Eigen INTERFACE -DEIGEN_USE_LAPACKE_STRICT)
    endif()

    target_compile_definitions(Eigen3::Eigen INTERFACE -DEIGEN_DONT_ALIGN_STATICALLY)
#    target_compile_definitions(Eigen3::Eigen INTERFACE -DEIGEN_DONT_VECTORIZE)
endif()


function(find_ceres)
    if(DMRG_PACKAGE_MANAGER MATCHES "find|cmake")
        foreach (tgt glog::glog;gflags;Eigen3::Eigen)
            if(NOT TARGET ${tgt})
                list(APPEND CERES_MISSING_TARGET ${tgt})
                mark_as_advanced(CERES_MISSING_TARGETS)
            endif()
        endforeach()
        if(CERES_MISSING_TARGET)
            message(FATAL_ERROR "Ceres: dependencies missing [${CERES_MISSING_TARGET}]")
        endif()
    endif()
    find_package(Ceres 2.0 PATH_SUFFIXES ceres ceres/lib ${N5} ${N6} ${N7} ${N8} ${REQUIRED}) # Flags ignore system packages. See cmake/SetupPaths.cmake)

    if(NOT Ceres_FOUND AND DMRG_PACKAGE_MANAGER MATCHES "cmake")
        message(STATUS "Ceres will be installed into ${DMRG_DEPS_INSTALL_DIR} on first build.")
        list(APPEND CERES_CMAKE_OPTIONS  -DEigen3_ROOT:PATH=${DMRG_DEPS_INSTALL_DIR})
        list(APPEND CERES_CMAKE_OPTIONS  -Dgflags_ROOT:PATH=${DMRG_DEPS_INSTALL_DIR})
        list(APPEND CERES_CMAKE_OPTIONS  -Dglog_ROOT:PATH=${DMRG_DEPS_INSTALL_DIR})
        message(STATUS "ceres options: ${CERES_CMAKE_OPTIONS}")
        include(cmake/InstallPackage.cmake)
        install_package(ceres "${DMRG_DEPS_INSTALL_DIR}" "${CERES_CMAKE_OPTIONS}" )
        find_package(Ceres HINTS ${DMRG_DEPS_INSTALL_DIR} NO_DEFAULT_PATH REQUIRED)
    endif()
    if(Ceres_FOUND AND TARGET Ceres::ceres)
        # Use this for the ceres targets defined by CONFIG mode find_package
        # These find_packages have the tendency to do the wrong thing, like
        #   - injecting shared libraries into static builds
        #   - using "-lpthread" instead of "pthread"
        get_target_property(CERES_TYPE ceres TYPE)
        if(CERES_TYPE MATCHES "SHARED" AND NOT BUILD_SHARED_LIBS)
            include(cmake/PrintTargetProperties.cmake)
            print_target_properties(Ceres::ceres)
            message(FATAL_ERROR "Target Ceres::ceres contains a shared library on a static build!")
        endif()
        #Remove any shared libraries like unwind etc which pollute static builds
        # Or... just relink it entirely
        include(cmake/TargetFilters.cmake)
        remove_library_shallow(Ceres::ceres "gcc_eh|unwind|lzma|Threads::Threads|pthread")
        target_link_libraries(Ceres::ceres INTERFACE pthread )

        # Check that ceres actually works
        include(cmake/CheckCeresCompiles.cmake)
        check_ceres_compiles("Ceres::ceres" "" "" "" "")
        if(NOT CERES_COMPILES)
            include(cmake/PrintTargetProperties.cmake)
            print_target_properties(Ceres::ceres)
            message(FATAL_ERROR "Could not compile simple ceres program")
        endif()

    else()
        message(FATAL_ERROR "Ceres could not be found")
    endif()
endfunction()


find_ceres()
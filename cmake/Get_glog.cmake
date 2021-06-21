
function(find_glog)
    if(DMRG_PACKAGE_MANAGER MATCHES "find|cmake")
        if(NOT TARGET gflags)
            message(FATAL_ERROR "glog::glog: dependencies missing [gflags]")
        endif()
    endif()
    if(DMRG_PACKAGE_MANAGER STREQUAL "find")
        set(REQUIRED REQUIRED)
    endif()
    if(DMRG_PACKAGE_MANAGER STREQUAL "cmake")
        set(NO_DEFAULT_PATH NO_DEFAULT_PATH)
    endif()
    # Glog can be compiled with or without libunwind
    # unwind can sometimes need gcc_eh and lzma as dependencies.
    # Yes, it's a mess. Here we add targets if they are present in the system
    # If not, we hope for the best making dummy targets.

    # ... Additionally, if only lzma is linked we get this weird error when linking ceres:
    #       undefined reference to
    #       `ceres::GradientProblemSolver::Summary::Summary()'
    # This could also be caused by using openblas from anaconda on Clang builds.
    #   UPDATE: The statement above may just be due to weird ABI incompatibility with conda libs sometimes.
    #           A better solution may be to avoid conda when that happens.

    find_package(Unwind) # If found defines target unwind::unwind
    find_library(LZMA_LIB NAMES lzma)
    find_library(GCC_EH_LIB NAMES gcc_eh)
    add_library(unwind::full INTERFACE IMPORTED)
    if(TARGET unwind::unwind AND LZMA_LIB AND GCC_EH_LIB)
        target_link_libraries(unwind::full INTERFACE gcc_eh unwind::unwind lzma)
    endif()

    find_package(glog 0.4
            HINTS ${glog_ROOT} ${DMRG_DEPS_INSTALL_DIR}
            ${NO_DEFAULT_PATH}
            ${REQUIRED})


    if(NOT glog_FOUND AND DMRG_PACKAGE_MANAGER MATCHES "cmake")
        message(STATUS "glog will be installed into ${DMRG_DEPS_INSTALL_DIR}")
        include(cmake/InstallPackage.cmake)
        list(APPEND GLOG_CMAKE_OPTIONS -Dgflags_ROOT:PATH=${DMRG_DEPS_INSTALL_DIR})
        install_package(glog "${DMRG_DEPS_INSTALL_DIR}" "${GLOG_CMAKE_OPTIONS}")
        find_package(glog 0.4
                HINTS ${glog_ROOT} ${DMRG_DEPS_INSTALL_DIR}
                NO_DEFAULT_PATH
                REQUIRED)
    endif()



    if(glog_FOUND AND TARGET glog::glog)
        get_target_property(GLOG_TYPE glog::glog TYPE)
        if(GLOG_TYPE MATCHES "SHARED" AND NOT BUILD_SHARED_LIBS)
            include(cmake/PrintTargetProperties.cmake)
            print_target_properties(glog::glog)
            message(FATAL_ERROR "Found shared glog library on a static build!")
        endif()

        include(cmake/TargetFilters.cmake)
        remove_library_shallow(glog::glog "Threads::Threads|pthread|unwind")
        target_link_libraries(glog::glog INTERFACE unwind::full pthread )
        include(cmake/CheckGlogCompiles.cmake)
        check_glog_compiles("glog::glog" "" "" "" "")
        if(NOT GLOG_COMPILES)
            include(cmake/PrintTargetProperties.cmake)
            print_target_properties(glog::glog)
            message(FATAL_ERROR "Could not compile a simple glog program")
        endif()
    else()
        message(FATAL_ERROR "glog not found")
    endif()
endfunction()

find_glog()
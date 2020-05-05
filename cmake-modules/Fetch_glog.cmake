if(DMRG_DOWNLOAD_METHOD MATCHES "find|fetch")
    if(NOT TARGET gflags::gflags)
        message(FATAL_ERROR "glog::glog: dependencies missing [gflags::gflags]")
    endif()
endif()




if(NOT TARGET glog::glog AND DMRG_DOWNLOAD_METHOD MATCHES "find|fetch")

    include(cmake-modules/CheckGlogCompiles.cmake)

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
    message(STATUS "Looking for glog config")
    # Glog should only look in conda on shared builds! Conda does not give us the static version
    if(DMRG_PREFER_CONDA_LIBS AND NOT BUILD_SHARED_LIBS)
        # Static case - conda should not be searched
        message(STATUS "Excluding conda from glog search in static builds")
        find_package(glog 0.4 HINTS ${CMAKE_INSTALL_PREFIX} PATHS $ENV{EBROOTGLOG} NO_CMAKE_PATH NO_CMAKE_ENVIRONMENT_PATH NO_CMAKE_PACKAGE_REGISTRY)
    else()
        find_package(glog 0.4 NO_CMAKE_PACKAGE_REGISTRY)
    endif()
    if(TARGET glog::glog)
        message(STATUS "Looking for glog config - found")
    else()
        message(STATUS "Looking for glog config - not found")
    endif()
    if(NOT TARGET glog::glog)
        message(STATUS "Looking for glog in system")
        find_library(GLOG_LIBRARIES     glog           )
        find_path   (GLOG_INCLUDE_DIR   glog/logging.h )
        if(GLOG_LIBRARIES AND GLOG_INCLUDE_DIR)
            check_glog_compiles("" "${GLOG_LIBRARIES};unwind::full;gflags::gflags;pthread" "${GLOG_INCLUDE_DIR}" "" "")
            if (GLOG_COMPILES)
                add_library(glog::glog ${LINK_TYPE} IMPORTED)
                set_target_properties(glog::glog PROPERTIES IMPORTED_LOCATION "${GLOG_LIBRARIES}")
                target_include_directories(glog::glog SYSTEM INTERFACE ${GLOG_INCLUDE_DIR})
                message(STATUS "Looking for glog in system - found: ${GLOG_LIBRARIES}")
            endif()
        else()
            message(STATUS "Looking for glog in system - not found")
        endif()

    endif()
endif()


if(NOT TARGET glog::glog AND DMRG_DOWNLOAD_METHOD MATCHES "fetch")
    message(STATUS "glog will be installed into ${CMAKE_INSTALL_PREFIX}")
    include(${PROJECT_SOURCE_DIR}/cmake-modules/BuildDependency.cmake)
    list(APPEND GLOG_CMAKE_OPTIONS -Dgflags_DIR:PATH=${CMAKE_INSTALL_PREFIX}/gflags/lib/cmake/gflags)
    list(APPEND GLOG_CMAKE_OPTIONS -DCMAKE_PREFIX_PATH:PATH=${CMAKE_INSTALL_PREFIX}/gflags)
    build_dependency(glog "${CMAKE_INSTALL_PREFIX}" "${GLOG_CMAKE_OPTIONS}")
    find_package(glog 0.4 NO_CMAKE_PACKAGE_REGISTRY)
    if(TARGET glog::glog)
        message(STATUS "glog installed successfully")
    else()
        message(FATAL_ERROR "glog could not be installed")
    endif()
endif()


if(TARGET glog::glog)
    get_target_property(GLOG_TYPE glog::glog TYPE)
    if(GLOG_TYPE MATCHES "SHARED" AND NOT BUILD_SHARED_LIBS)
        include(cmake-modules/PrintTargetProperties.cmake)
        print_target_properties(glog::glog)
        message(FATAL_ERROR "Found shared glog library on a static build!")
    endif()

    include(cmake-modules/TargetFilters.cmake)
    remove_library_shallow(glog::glog "Threads::Threads|pthread|unwind|gflags")
    target_link_libraries(glog::glog INTERFACE unwind::full gflags::gflags pthread )

    # Modernize
    get_property(imp_loc_set TARGET glog::glog PROPERTY IMPORTED_LOCATION SET) # Returns a boolean if set
    get_property(loc_set     TARGET glog::glog PROPERTY LOCATION SET) # Returns a boolean if set
    if(loc_set AND NOT imp_loc_set)
        get_target_property(imp_loc glog::glog LOCATION)
        set_target_properties(glog::glog PROPERTIES IMPORTED_LOCATION ${imp_loc})
    endif()

    check_glog_compiles("glog::glog" "" "" "${GLOG_FLAGS}" "")
    if(NOT GLOG_COMPILES)
        include(cmake-modules/PrintTargetProperties.cmake)
        print_target_properties(glog::glog)
        message(FATAL_ERROR "Could not compile a simple glog program")
    endif()
endif()
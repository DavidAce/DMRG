# Gflags comes in static flavor in conda also!
set(GFLAGS_HINTS $ENV{EBROOTGFLAGS} ${CMAKE_INSTALL_PREFIX} ${CONDA_HINTS})


find_package(gflags COMPONENTS nothreads_static HINTS ${GFLAGS_HINTS}
        PATHS $ENV{EBROOTGFLAGS}
        PATH_SUFFIXES gflags gflags/lib)


if(TARGET gflags)
    message(STATUS "gflags found")

elseif(DOWNLOAD_MISSING)
    message(STATUS "gflags will be installed into ${CMAKE_INSTALL_PREFIX}")
    include(${PROJECT_SOURCE_DIR}/cmake-modules/BuildDependency.cmake)
    build_dependency(gflags "${CMAKE_INSTALL_PREFIX}" "")
    find_package(gflags
            HINTS ${CMAKE_INSTALL_PREFIX}
            PATH_SUFFIXES gflags gflags/lib
            NO_DEFAULT_PATH)
    if(TARGET gflags)
        message(STATUS "gflags installed successfully")
    else()
        message(FATAL_ERROR "gflags could not be downloaded.")
    endif()

else()
    message(FATAL_ERROR "Dependency gflags not found and DOWNLOAD_MISSING is OFF")
endif()

if(TARGET gflags)
    include(cmake-modules/filterTarget.cmake)
    remove_shared(gflags)
    remove_pthread_shallow(gflags)

    get_target_property(GFLAGS_TYPE gflags TYPE)
    if(GFLAGS_TYPE MATCHES "SHARED" AND NOT BUILD_SHARED_LIBS)
        include(cmake-modules/PrintTargetProperties.cmake)
        print_target_properties(gflags)
        message(FATAL_ERROR "Found shared gflags library on a static build!")
    endif()

    if(NOT TARGET gflags::gflags)
        # Copy gflags to gflags::gflags to follow proper naming convention
        include(cmake-modules/CopyTarget.cmake)
        copy_target(gflags::gflags gflags)
    endif()

endif()
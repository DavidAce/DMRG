
if(NOT BUILD_SHARED_LIBS AND NOT DMRG_PREFER_CONDA_LIBS)
    set(COMPONENTS COMPONENTS)
    set(ITEMS nothreads_static)
endif()

if(BUILD_SHARED_LIBS AND DMRG_PREFER_CONDA_LIBS)
    list(APPEND GFLAGS_CONDA_PREFIX ${DMRG_CONDA_CANDIDATE_PATHS})
elseif(DMRG_PREFER_CONDA_LIBS)
    # Static build should avoid conda gflags since it only offers a shared version
    set(NO_DEFAULT_PATH NO_DEFAULT_PATH)
endif()

if(NOT TARGET gflags AND DMRG_PACKAGE_MANAGER MATCHES "find|cmake")
    # Gflags comes in static flavor in conda also!
    find_package(gflags
            HINTS ${GFLAGS_CONDA_PREFIX}
            ${COMPONENTS} ${ITEMS}
            NO_CMAKE_PACKAGE_REGISTRY
            ${NO_DEFAULT_PATH}
            )
    if(TARGET gflags)
        message(STATUS "Found gflags")
    endif()
endif()

if(NOT TARGET gflags AND DMRG_PACKAGE_MANAGER MATCHES "cmake" )
    message(STATUS "gflags will be installed into ${CMAKE_INSTALL_PREFIX}")
    include(${PROJECT_SOURCE_DIR}/cmake-modules/BuildDependency.cmake)
    build_dependency(gflags "${CMAKE_INSTALL_PREFIX}" "")
    find_package(gflags
            HINTS ${GFLAGS_CONDA_PREFIX}
            ${COMPONENTS} ${ITEMS}
            HINTS ${CMAKE_INSTALL_PREFIX}/gflags
            NO_CMAKE_PACKAGE_REGISTRY
            ${NO_DEFAULT_PATH}
            )
    if(TARGET gflags)
        message(STATUS "gflags installed successfully")
    else()
        message(FATAL_ERROR "gflags could not be installed")
    endif()
endif()

if(TARGET gflags)
    get_target_property(GFLAGS_TYPE gflags TYPE)
    if(GFLAGS_TYPE MATCHES "SHARED" AND NOT BUILD_SHARED_LIBS)
        include(cmake-modules/PrintTargetProperties.cmake)
        print_target_properties(gflags)
        message(STATUS "${CMAKE_PREFIX_PATH}" )
        message(FATAL_ERROR "Target gflags contains a shared library on a static build!")
    endif()

    if(NOT BUILD_SHARED_LIBS)
        include(cmake-modules/TargetFilters.cmake)
        replace_or_remove_shared(gflags)
        replace_pthread_shallow(gflags)
    endif()

    # Modernize
    get_property(imp_loc_set TARGET gflags PROPERTY IMPORTED_LOCATION SET) # Returns a boolean if set
    get_property(loc_set     TARGET gflags PROPERTY LOCATION SET) # Returns a boolean if set
    if(loc_set AND NOT imp_loc_set)
        get_target_property(imp_loc gflags LOCATION)
        set_target_properties(gflags PROPERTIES IMPORTED_LOCATION ${imp_loc})
    endif()

#    if(NOT TARGET gflags::gflags)
        # Copy gflags to gflags::gflags to follow proper naming convention
#        add_library(gflags::gflags ALIAS gflags)
#        include(cmake-modules/CopyTarget.cmake)
#        copy_target(gflags gflags::gflags)
#    endif()
endif()
function(find_gflags)
    if(DMRG_PACKAGE_MANAGER STREQUAL "find")
        set(REQUIRED REQUIRED)
    endif()
    if(DMRG_PACKAGE_MANAGER STREQUAL "cmake")
        set(NO_DEFAULT_PATH NO_DEFAULT_PATH)
    endif()
    if(NOT BUILD_SHARED_LIBS)
        set(COMPONENTS COMPONENTS)
        set(ITEMS nothreads_static)
    endif()

    # Gflags comes in static flavor in conda also!
    find_package(gflags
            HINTS ${gflags_ROOT} ${DMRG_DEPS_INSTALL_DIR}
            ${COMPONENTS} ${ITEMS}
            ${NO_DEFAULT_PATH}
            ${REQUIRED}
            )
    if(NOT gflags_FOUND AND DMRG_PACKAGE_MANAGER MATCHES "cmake" )
        message(STATUS "gflags will be installed into ${DMRG_DEPS_INSTALL_DIR}")
        include(${PROJECT_SOURCE_DIR}/cmake/InstallPackage.cmake)
        install_package(gflags "${DMRG_DEPS_INSTALL_DIR}" "")
        find_package(gflags
                HINTS ${GFLAGS_CONDA_PREFIX}
                ${COMPONENTS} ${ITEMS}
                HINTS ${gflags_ROOT} ${DMRG_DEPS_INSTALL_DIR}
                NO_DEFAULT_PATH
                REQUIRED)
    endif()
    if(gflags_FOUND AND TARGET gflags)
        get_target_property(GFLAGS_TYPE gflags TYPE)
        if(GFLAGS_TYPE MATCHES "SHARED" AND NOT BUILD_SHARED_LIBS)
            include(cmake/PrintTargetProperties.cmake)
            print_target_properties(gflags)
            message(STATUS "${CMAKE_PREFIX_PATH}" )
            message(FATAL_ERROR "Target gflags contains a shared library on a static build!")
        endif()

        if(NOT BUILD_SHARED_LIBS)
            include(cmake/TargetFilters.cmake)
            replace_or_remove_shared(gflags)
            replace_pthread_shallow(gflags)
        endif()

        # Modernize
#        get_property(imp_loc_set TARGET gflags PROPERTY IMPORTED_LOCATION SET) # Returns a boolean if set
#        get_property(loc_set     TARGET gflags PROPERTY LOCATION SET) # Returns a boolean if set
#        if(loc_set AND NOT imp_loc_set)
#            get_target_property(imp_loc gflags LOCATION)
#            set_target_properties(gflags PROPERTIES IMPORTED_LOCATION ${imp_loc})
#        endif()
    else()
        message(FATAL_ERROR "gflags not found")
    endif()

    if(NOT TARGET gflags::gflags)
        add_library(gflags::gflags ALIAS gflags)
    endif()

endfunction()
find_gflags()
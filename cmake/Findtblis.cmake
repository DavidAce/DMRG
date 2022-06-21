unset(TBLIS_LIBRARY)
unset(TBLIS_LIBRARY CACHE)
find_library(TBLIS_LIBRARY
        tblis
        HINTS ${DMRG_DEPS_INSTALL_DIR}
        PATH_SUFFIXES lib tblis/lib
        NO_CMAKE_ENVIRONMENT_PATH
        NO_SYSTEM_ENVIRONMENT_PATH
        NO_CMAKE_SYSTEM_PATH
        )
find_path(TBLIS_INCLUDE_DIR
        tblis/tblis.h
        HINTS ${DMRG_DEPS_INSTALL_DIR}
        PATH_SUFFIXES include tblis/include
        NO_CMAKE_ENVIRONMENT_PATH
        NO_SYSTEM_ENVIRONMENT_PATH
        NO_CMAKE_SYSTEM_PATH
        )



include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(tblis
        DEFAULT_MSG
        TBLIS_LIBRARY TBLIS_INCLUDE_DIR)


if (tblis_FOUND)
    add_library(tblis::tblis UNKNOWN IMPORTED)
    set_target_properties(tblis::tblis PROPERTIES IMPORTED_LOCATION "${TBLIS_LIBRARY}")
    target_include_directories(tblis::tblis SYSTEM INTERFACE ${TBLIS_INCLUDE_DIR})
endif()
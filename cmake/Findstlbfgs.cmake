unset(STLBFGS_LIBRARY)
unset(STLBFGS_LIBRARY CACHE)
unset(STLBFGS_INCLUDE_DIR)
unset(STLBFGS_INCLUDE_DIR CACHE)
unset(STLBFGS_LS_INCLUDE_DIR)
unset(STLBFGS_LS_INCLUDE_DIR CACHE)

find_library(STLBFGS_LIBRARY
        stlbfgs
        HINTS ${DMRG_DEPS_INSTALL_DIR}
        PATH_SUFFIXES lib stlbfgs/lib
        NO_CMAKE_ENVIRONMENT_PATH
        NO_SYSTEM_ENVIRONMENT_PATH
        NO_CMAKE_SYSTEM_PATH
        )
find_path(STLBFGS_INCLUDE_DIR
        stlbfgs/stlbfgs.h
        HINTS ${DMRG_DEPS_INSTALL_DIR}
        PATH_SUFFIXES include stlbfgs/include
        NO_CMAKE_ENVIRONMENT_PATH
        NO_SYSTEM_ENVIRONMENT_PATH
        NO_CMAKE_SYSTEM_PATH
        )
find_path(STLBFGS_LS_INCLUDE_DIR
        stlbfgs/linesearch.h
        HINTS ${DMRG_DEPS_INSTALL_DIR}
        PATH_SUFFIXES include stlbfgs/include
        NO_CMAKE_ENVIRONMENT_PATH
        NO_SYSTEM_ENVIRONMENT_PATH
        NO_CMAKE_SYSTEM_PATH
        )


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(stlbfgs
        DEFAULT_MSG
        STLBFGS_LIBRARY STLBFGS_INCLUDE_DIR STLBFGS_LS_INCLUDE_DIR)


if (stlbfgs_FOUND)
    add_library(stlbfgs::stlbfgs UNKNOWN IMPORTED)
    set_target_properties(stlbfgs::stlbfgs PROPERTIES IMPORTED_LOCATION "${STLBFGS_LIBRARY}")
    target_include_directories(stlbfgs::stlbfgs SYSTEM INTERFACE ${STLBFGS_INCLUDE_DIR})
    target_compile_definitions(stlbfgs::stlbfgs INTERFACE STLBFGS_INT_SIZE=32)
endif ()
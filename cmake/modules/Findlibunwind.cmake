function(find_unwind)
    find_library(LIBUNWIND_LIBRARY
            unwind
            HINTS ${DMRG_DEPS_INSTALL_DIR}
            PATH_SUFFIXES lib libunwind/lib
            NO_CMAKE_ENVIRONMENT_PATH
            NO_SYSTEM_ENVIRONMENT_PATH
            NO_CMAKE_SYSTEM_PATH
            )
    find_path(LIBUNWIND_INCLUDE_DIR
            libunwind.h
            HINTS ${DMRG_DEPS_INSTALL_DIR}
            PATH_SUFFIXES include libunwind/include
            NO_CMAKE_ENVIRONMENT_PATH
            NO_SYSTEM_ENVIRONMENT_PATH
            NO_CMAKE_SYSTEM_PATH
            )

endfunction()

find_unwind()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(libunwind
        DEFAULT_MSG
        LIBUNWIND_LIBRARY LIBUNWIND_INCLUDE_DIR)


if (libunwind_FOUND AND NOT TARGET libunwind::libunwind)
    add_library(libunwind::libunwind STATIC IMPORTED)
    set_target_properties(libunwind::libunwind PROPERTIES IMPORTED_LOCATION "${LIBUNWIND_LIBRARY}")
    target_include_directories(libunwind::libunwind SYSTEM INTERFACE ${LIBUNWIND_INCLUDE_DIR})
endif ()
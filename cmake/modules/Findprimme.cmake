function(find_primme)
    if(NOT BUILD_SHARED_LIBS)
        set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX} ${CMAKE_SHARED_LIBRARY_SUFFIX})
    endif()
    unset(PRIMME_LIBRARY)
    unset(PRIMME_LIBRARY CACHE)
    find_library(PRIMME_LIBRARY
                 primme
                 HINTS ${DMRG_DEPS_INSTALL_DIR}
                 PATH_SUFFIXES lib primme/lib
                 NO_CMAKE_ENVIRONMENT_PATH
                 NO_SYSTEM_ENVIRONMENT_PATH
                 NO_CMAKE_SYSTEM_PATH
                 )
    find_path(PRIMME_INCLUDE_DIR
              primme/primme.h
              HINTS ${DMRG_DEPS_INSTALL_DIR}
              PATH_SUFFIXES include primme/include
              NO_CMAKE_ENVIRONMENT_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH
              )

endfunction()
find_primme()
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(primme
                                  DEFAULT_MSG
                                  PRIMME_LIBRARY PRIMME_INCLUDE_DIR)

if(primme_FOUND AND NOT TARGET primme::primme)
    add_library(primme::primme UNKNOWN IMPORTED)
    set_target_properties(primme::primme PROPERTIES IMPORTED_LOCATION "${PRIMME_LIBRARY}")
    target_include_directories(primme::primme SYSTEM INTERFACE ${PRIMME_INCLUDE_DIR})
    target_compile_definitions(primme::primme INTERFACE PRIMME_INT_SIZE=32)
endif()
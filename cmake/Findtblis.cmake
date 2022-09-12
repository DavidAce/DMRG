
function(find_tblis)
    find_library(TBLIS_LIBRARY
                 tblis
                 HINTS ${DMRG_DEPS_INSTALL_DIR} ${CMAKE_INSTALL_PREFIX}
                 PATH_SUFFIXES lib tblis/lib
                 NO_CMAKE_ENVIRONMENT_PATH
                 NO_SYSTEM_ENVIRONMENT_PATH
                 NO_CMAKE_SYSTEM_PATH
                 QUIET
                 )
    find_path(TBLIS_INCLUDE_DIR
              tblis/tblis.h
              HINTS ${DMRG_DEPS_INSTALL_DIR} ${CMAKE_INSTALL_PREFIX}
              PATH_SUFFIXES include tblis/include
              NO_CMAKE_ENVIRONMENT_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH
              QUIET
              )
    find_library(TCI_LIBRARY
                 tci
                 HINTS ${DMRG_DEPS_INSTALL_DIR} ${CMAKE_INSTALL_PREFIX}
                 PATH_SUFFIXES lib tblis/lib
                 NO_CMAKE_ENVIRONMENT_PATH
                 NO_SYSTEM_ENVIRONMENT_PATH
                 NO_CMAKE_SYSTEM_PATH
                 QUIET
                 )
    find_path(TCI_INCLUDE_DIR
              tci/tci_global.h
              HINTS ${DMRG_DEPS_INSTALL_DIR} ${CMAKE_INSTALL_PREFIX}
              PATH_SUFFIXES include tci/include
              NO_CMAKE_ENVIRONMENT_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH
              QUIET
              )

endfunction()

find_tblis()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(tblis
                                  DEFAULT_MSG
                                  TBLIS_LIBRARY TBLIS_INCLUDE_DIR TCI_LIBRARY TCI_INCLUDE_DIR)

if(tblis_FOUND AND NOT TARGET tblis::tblis)
    add_library(tblis::tblis UNKNOWN IMPORTED)
    add_library(tblis::tci UNKNOWN IMPORTED)
    set_target_properties(tblis::tci PROPERTIES IMPORTED_LOCATION "${TCI_LIBRARY}")
    set_target_properties(tblis::tblis PROPERTIES IMPORTED_LOCATION "${TBLIS_LIBRARY}")
    target_include_directories(tblis::tci SYSTEM INTERFACE ${TCI_INCLUDE_DIR})
    target_include_directories(tblis::tblis SYSTEM INTERFACE ${TBLIS_INCLUDE_DIR})
    target_link_libraries(tblis::tblis INTERFACE tblis::tci)
    if(NOT BUILD_SHARED_LIBS)
        target_link_libraries(tblis::tblis INTERFACE atomic hwloc)
    endif()
endif()
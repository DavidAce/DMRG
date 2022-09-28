
function(find_tblis)
    if(NOT BUILD_SHARED_LIBS)
        set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX} ${CMAKE_SHARED_LIBRARY_SUFFIX})
    endif()
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
function(check_tblis_compiles)
    set(CMAKE_REQUIRED_LIBRARIES tblis::tblis)
    if(NOT BUILD_SHARED_LIBS OR CMAKE_LINK_SEARCH_START_STATIC)
        set(CMAKE_REQUIRED_LINK_OPTIONS -static-libgcc)
    endif()
    check_cxx_source_compiles("
                        #include <tblis/tblis.h>
                        #include <tblis/util/thread.h>
                        int main() {
                           tblis::len_vector da;
                           return 0;
                        }
                        " COMPILES_TBLIS)
    if(NOT COMPILES_TBLIS)
        message(FATAL_ERROR "Failed to compile a simple program using tblis::tblis")
    endif()
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
    find_package(OpenMP COMPONENTS CXX REQUIRED)
    target_link_libraries(tblis::tci INTERFACE OpenMP::OpenMP_CXX)
    if(NOT BUILD_SHARED_LIBS)
        find_library(ATOMIC_LIBRARY atomic)
        find_library(HWLOC_LIBRARY hwloc)
        find_library(UDEV_LIBRARY udev)
        if(ATOMIC_LIBRARY)
            target_link_libraries(tblis::tblis INTERFACE ${ATOMIC_LIBRARY})
        else()
            target_link_libraries(tblis::tblis INTERFACE atomic)
            message(WARNING "atomic library not found (tblis may depend on it). Linking -latomic")
        endif()
        if(HWLOC_LIBRARY)
            target_link_libraries(tblis::tblis INTERFACE ${HWLOC_LIBRARY})
        else()
            message(WARNING "hwloc library not found (tblis may depend on it). Linking -lhwloc")
            target_link_libraries(tblis::tblis INTERFACE hwloc)
        endif()
        if(UDEV_LIBRARY)
            target_link_libraries(tblis::tblis INTERFACE ${UDEV_LIBRARY})
        else()
            message(WARNING "udev library not found (tblis may depend on it). Linking -ludev")
            target_link_libraries(tblis::tblis INTERFACE udev)
        endif()
    endif()
    target_link_libraries(tblis::tblis INTERFACE tblis::tci)

    check_tblis_compiles()
endif()


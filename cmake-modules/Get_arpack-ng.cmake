if(DMRG_PACKAGE_MANAGER MATCHES "find|cmake")
    foreach (tgt BLAS::BLAS;LAPACK::LAPACK;gfortran::gfortran)
        if(NOT TARGET ${tgt})
            list(APPEND arpack_ng_MISSING_TARGET ${tgt})
            mark_as_advanced(arpack_ng_MISSING_TARGET)
        endif()
    endforeach()
    if(arpack_ng_MISSING_TARGET)
        message(FATAL_ERROR "arpack-ng: dependencies missing [${arpack_ng_MISSING_TARGET}]")
    endif()
endif()



mark_as_advanced(ARPACK_DIRECTORY_HINTS)

function(find_arpack_ng REQUIRED)
    unset(arpack_ng_LIBRARIES)
    unset(arpack_ng_INCLUDE_DIRS)
    unset(arpack_ng_LIBRARIES CACHE)
    unset(arpack_ng_INCLUDE_DIRS CACHE)

    list(APPEND ARPACK_DIRECTORY_HINTS
            ${CMAKE_INSTALL_PREFIX}/arpack-ng
            ${CMAKE_INSTALL_PREFIX}/arpack-ng/lib
            ${CMAKE_INSTALL_PREFIX}/arpack-ng/lib64
            ${CMAKE_INSTALL_PREFIX}
            ${CMAKE_INSTALL_PREFIX}/lib
            ${CMAKE_INSTALL_PREFIX}/lib64
            )
    find_library(arpack_ng_LIBRARIES
            arpack
            HINTS ${ARPACK_DIRECTORY_HINTS}
            ${REQUIRED}
            NO_CMAKE_PACKAGE_REGISTRY
            )
    find_path(arpack_ng_INCLUDE_DIRS
            arpack/debug.h
            HINTS ${ARPACK_DIRECTORY_HINTS}
            PATH_SUFFIXES include
            NO_CMAKE_PACKAGE_REGISTRY
            )
    if(arpack_ng_LIBRARIES)
        set(arpack_ng_LIBRARIES PARENT_SCOPE)
    endif()
    if(arpack_ng_INCLUDE_DIRS)
        set(arpack_ng_INCLUDE_DIRS PARENT_SCOPE)
    endif()
endfunction()





if(NOT TARGET ARPACK::ARPACK AND DMRG_PACKAGE_MANAGER MATCHES "find|cmake")
    find_package(arpack-ng QUIET)
    if(TARGET ARPACK::ARPACK)
        message(STATUS "Successfully installed arpack-ng")
        target_link_libraries(ARPACK::ARPACK INTERFACE gfortran::gfortran)
        target_link_libraries(ARPACK::ARPACK INTERFACE BLAS::BLAS LAPACK::LAPACK gfortran::gfortran)
    endif()
endif()

if (NOT TARGET ARPACK::ARPACK AND DMRG_PACKAGE_MANAGER STREQUAL "find")
    message(STATUS "Looking for arpack-ng in system")
    find_arpack_ng(REQUIRED)
    if(NOT arpack_ng_LIBRARIES)
        message(STATUS "Looking for arpack-ng in system - not found")
    else()
        message(STATUS "Looking for arpack-ng in system - found: ${arpack_ng_LIBRARIES}")
        add_library(ARPACK::ARPACK ${LINK_TYPE} IMPORTED)
        set_target_properties(ARPACK::ARPACK PROPERTIES IMPORTED_LOCATION "${arpack_ng_LIBRARIES}")
        target_link_libraries(ARPACK::ARPACK INTERFACE BLAS::BLAS LAPACK::LAPACK gfortran::gfortran)
        if("${arpack_ng_LIBRARIES}" MATCHES "usr" AND BUILD_SHARED_LIBS)
            target_link_libraries(ARPACK::ARPACK INTERFACE BLAS::BLAS )
            message(WARNING "Found arpack-ng in system as shared library. "
                            "Make sure to use correct libraries libblas.so.3 and liblapack.so.3 "
                            "by using \n"
                            "sudo update-alternatives --config libblas.so.3-x86_64-linux-gnu \n"
                            "sudo update-alternatives --config liblapack.so.3-x86_64-linux-gnu ")
        endif()
    endif()
endif()


if(NOT TARGET ARPACK::ARPACK AND DMRG_PACKAGE_MANAGER MATCHES "cmake")
    message(STATUS "arpack-ng will be installed into ${CMAKE_INSTALL_PREFIX}/arpack-ng")
    include(${PROJECT_SOURCE_DIR}/cmake-modules/BuildDependency.cmake)
    build_dependency(arpack-ng "${CMAKE_INSTALL_PREFIX}" "${ARPACK_CMAKE_OPTIONS}")
    find_package(arpack-ng HINTS ${CMAKE_INSTALL_PREFIX}/arpack-ng REQUIRED NO_DEFAULT_PATH)
    if(TARGET ARPACK::ARPACK)
        message(STATUS "Successfully installed arpack-ng")
        target_link_libraries(ARPACK::ARPACK INTERFACE BLAS::BLAS LAPACK::LAPACK gfortran::gfortran)
    else()
        message(FATAL_ERROR "arpack-ng could not be installed")
    endif()
endif()
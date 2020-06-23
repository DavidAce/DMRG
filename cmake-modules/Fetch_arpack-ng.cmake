if(DMRG_DOWNLOAD_METHOD MATCHES "find|fetch")
    foreach (tgt blas::blas;lapack::lapack;gfortran::gfortran)
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





if(NOT TARGET arpack::arpack AND DMRG_DOWNLOAD_METHOD MATCHES "find|fetch")
    find_arpack_ng("")
    if(arpack_ng_LIBRARIES AND arpack_ng_INCLUDE_DIRS)
        message(STATUS "Found arpack-ng")
        add_library(arpack::arpack ${LINK_TYPE} IMPORTED)
        set_target_properties(arpack::arpack PROPERTIES IMPORTED_LOCATION "${arpack_ng_LIBRARIES}")
        target_include_directories(arpack::arpack SYSTEM INTERFACE ${arpack_ng_INCLUDE_DIRS})
        target_link_libraries(arpack::arpack INTERFACE blas::blas lapack::lapack gfortran::gfortran)
    endif()
endif()

if (NOT TARGET arpack::arpack AND DMRG_DOWNLOAD_METHOD STREQUAL "find")
    message(STATUS "Looking for arpack-ng in system")
    find_arpack_ng(REQUIRED)
    if(NOT arpack_ng_LIBRARIES)
        message(STATUS "Looking for arpack-ng in system - not found")
    else()
        message(STATUS "Looking for arpack-ng in system - found: ${arpack_ng_LIBRARIES}")
        add_library(arpack::arpack ${LINK_TYPE} IMPORTED)
        set_target_properties(arpack::arpack PROPERTIES IMPORTED_LOCATION "${arpack_ng_LIBRARIES}")
        target_link_libraries(arpack::arpack INTERFACE blas::blas lapack::lapack gfortran::gfortran)
        if("${arpack_ng_LIBRARIES}" MATCHES "usr" AND BUILD_SHARED_LIBS)
            target_link_libraries(arpack::arpack INTERFACE blas::blas )
            message(WARNING "Found arpack-ng in system as shared library. "
                            "Make sure to use correct libraries libblas.so.3 and liblapack.so.3 "
                            "by using \n"
                            "sudo update-alternatives --config libblas.so.3-x86_64-linux-gnu \n"
                            "sudo update-alternatives --config liblapack.so.3-x86_64-linux-gnu ")
        endif()
    endif()
endif()


if(NOT TARGET arpack::arpack AND DMRG_DOWNLOAD_METHOD MATCHES "fetch")
    message(STATUS "arpack-ng will be installed into ${CMAKE_INSTALL_PREFIX}/arpack-ng")
    #####################################################################
    ### Prepare lists with generator expressions, replacing all semicolons.
    ### Otherwise, passing raw lists results in only the first element
    ### of the list to be passed.
    ####################################################################
    include(${PROJECT_SOURCE_DIR}/cmake-modules/getExpandedTarget.cmake)
    expand_target_libs(blas::blas       EXPANDED_BLAS)
    expand_target_libs(lapack::lapack   EXPANDED_LAPACK)
    string (REPLACE ";" "$<SEMICOLON>" EXPANDED_BLAS_GENERATOR      "${EXPANDED_BLAS}")
    string (REPLACE ";" "$<SEMICOLON>" EXPANDED_LAPACK_GENERATOR    "${EXPANDED_LAPACK}")
    set(ARPACK_FLAGS "-w -m64 -fPIC")
    include(${PROJECT_SOURCE_DIR}/cmake-modules/BuildDependency.cmake)
    list(APPEND ARPACK_CMAKE_OPTIONS   -DARPACK_FLAGS=${ARPACK_FLAGS})
    list(APPEND ARPACK_CMAKE_OPTIONS   -DBLAS_LIBRARIES=${EXPANDED_BLAS_GENERATOR})
    list(APPEND ARPACK_CMAKE_OPTIONS   -DLAPACK_LIBRARIES=${EXPANDED_LAPACK_GENERATOR})
    list(APPEND ARPACK_CMAKE_OPTIONS   -DCMAKE_VERBOSE_MAKEFILE=OFF)
    mark_as_advanced(ARPACK_CMAKE_OPTIONS)
    mark_as_advanced(ARPACK_FLAGS)
    build_dependency(arpack-ng "${CMAKE_INSTALL_PREFIX}" "${ARPACK_CMAKE_OPTIONS}")
    find_arpack_ng(REQUIRED)
    if(arpack_ng_LIBRARIES AND arpack_ng_INCLUDE_DIRS)
        message(STATUS "Successfully installed arpack-ng")
        add_library(arpack::arpack ${LINK_TYPE} IMPORTED)
        set_target_properties(arpack::arpack PROPERTIES IMPORTED_LOCATION "${arpack_ng_LIBRARIES}")
        target_include_directories(arpack::arpack SYSTEM INTERFACE ${arpack_ng_INCLUDE_DIRS})
        target_link_libraries(arpack::arpack INTERFACE blas::blas lapack::lapack gfortran::gfortran)
        mark_as_advanced(arpack_ng_LIBRARIES)
        mark_as_advanced(arpack_ng_INCLUDE_DIRS)
    else()
        message(FATAL_ERROR "arpack-ng could not be installed")
    endif()
endif()
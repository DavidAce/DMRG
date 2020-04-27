if(NOT TARGET arpack::arpack AND DMRG_DOWNLOAD_METHOD MATCHES "find|fetch|native|conan")
    include(GNUInstallDirs)
    unset(arpack_ng_LIBRARIES)
    unset(arpack_ng_INCLUDE_DIRS)
    find_package(arpack-ng HINTS ${CMAKE_INSTALL_PREFIX} PATH_SUFFIXES lib lib/cmake arpack-ng/${CMAKE_INSTALL_LIBDIR} NO_DEFAULT_PATH )
    if(arpack_ng_LIBRARIES AND arpack_ng_INCLUDE_DIRS)
        message(STATUS "Found arpack-ng")
        add_library(arpack::arpack ${LINK_TYPE} IMPORTED)
        set_target_properties(arpack::arpack PROPERTIES IMPORTED_LOCATION "${arpack_ng_LIBRARIES}")
        target_include_directories(arpack::arpack SYSTEM INTERFACE ${arpack_ng_INCLUDE_DIRS})
        target_link_libraries(arpack::arpack INTERFACE blas::blas lapack::lapack gfortran::gfortran)
    endif()

    if (NOT TARGET arpack::arpack AND DMRG_DOWNLOAD_METHOD MATCHES "find")
        message(STATUS "Searching for arpack-ng in system")
        find_library(ARPACK_LIBRARIES
                NAMES arpack
                HINTS ${CMAKE_INSTALL_PREFIX} $ENV{EBROOTARPACKMINNG}
                PATH_SUFFIXES arpack-ng arpack arpack/lib arpack-ng/lib
                )
        if(NOT ARPACK_LIBRARIES)
            message(STATUS "Searching for arpack-ng - failed")
        else()
            message(STATUS "Searching for arpack-ng - Success: ${ARPACK_LIBRARIES}")
            add_library(arpack::arpack ${LINK_TYPE} IMPORTED)
            set_target_properties(arpack::arpack PROPERTIES IMPORTED_LOCATION "${ARPACK_LIBRARIES}")
            target_link_libraries(arpack::arpack INTERFACE blas::blas lapack::lapack  gfortran::gfortran)
            if("${ARPACK_LIBRARIES}" MATCHES "usr" AND BUILD_SHARED_LIBS)
                target_link_libraries(arpack::arpack INTERFACE openblas::openblas )
                message(WARNING "Found arpack-ng in system as shared library. "
                                "Make sure to use correct libraries libblas.so.3 and liblapack.so.3 "
                                "by using \n"
                                "sudo update-alternatives --config libblas.so.3-x86_64-linux-gnu \n"
                                "sudo update-alternatives --config liblapack.so.3-x86_64-linux-gnu ")
            endif()
        endif()
    endif()
endif()


if(NOT TARGET arpack::arpack AND DMRG_DOWNLOAD_METHOD MATCHES "fetch|native|conan")
    message(STATUS "Arpack-ng will be installed into ${CMAKE_INSTALL_PREFIX}/arpack-ng")
    if(BUILD_SHARED_LIBS)
        set(ARPACK_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
    else()
        set(ARPACK_SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})
    endif()

    #####################################################################
    ### Prepare lists with generator expressions, replacing all semicolons.
    ### Otherwise, passing raw lists results  in only the first element
    ### of the list to be passed.
    ####################################################################
    include(${PROJECT_SOURCE_DIR}/cmake-modules/getExpandedTarget.cmake)
    include(${PROJECT_SOURCE_DIR}/cmake-modules/TargetFilters.cmake)
    include(${PROJECT_SOURCE_DIR}/cmake-modules/PrintTargetInfo.cmake)
    expand_target_libs(blas::blas       EXPANDED_BLAS)
    expand_target_libs(lapack::lapack   EXPANDED_LAPACK)
    string (REPLACE ";" "$<SEMICOLON>" EXPANDED_BLAS_GENERATOR      "${EXPANDED_BLAS}")
    string (REPLACE ";" "$<SEMICOLON>" EXPANDED_LAPACK_GENERATOR    "${EXPANDED_LAPACK}")

    ####################################################################
    set(ARPACK_FLAGS "-w -m64 -fPIC")
    message(STATUS "arpack-ng will be installed into ${CMAKE_INSTALL_PREFIX}")
    include(${PROJECT_SOURCE_DIR}/cmake-modules/BuildDependency.cmake)
    list(APPEND ARPACK_CMAKE_OPTIONS   -DARPACK_FLAGS=${ARPACK_FLAGS})
    list(APPEND ARPACK_CMAKE_OPTIONS   -DBLAS_LIBRARIES=${EXPANDED_BLAS_GENERATOR})
    list(APPEND ARPACK_CMAKE_OPTIONS   -DLAPACK_LIBRARIES=${EXPANDED_LAPACK_GENERATOR})
    mark_as_advanced(ARPACK_CMAKE_OPTIONS)
    build_dependency(arpack-ng "${CMAKE_INSTALL_PREFIX}/arpack-ng" "${ARPACK_CMAKE_OPTIONS}")

    find_package(arpack-ng HINTS ${CMAKE_INSTALL_PREFIX} PATH_SUFFIXES lib lib/cmake arpack-ng/${CMAKE_INSTALL_LIBDIR} REQUIRED NO_DEFAULT_PATH )
    if(arpack_ng_LIBRARIES AND arpack_ng_INCLUDE_DIRS)
        message(STATUS "Successfully installed arpack-ng")
        message(STATUS "arpack_ng_LIBRARIES     : ${arpack_ng_LIBRARIES}")
        message(STATUS "arpack_ng_INCLUDE_DIRS  : ${arpack_ng_INCLUDE_DIRS}")
        add_library(arpack::arpack ${LINK_TYPE} IMPORTED)
        set_target_properties(arpack::arpack PROPERTIES IMPORTED_LOCATION "${arpack_ng_LIBRARIES}")
        target_include_directories(arpack::arpack SYSTEM INTERFACE ${arpack_ng_INCLUDE_DIRS})
        target_link_libraries(arpack::arpack INTERFACE blas::blas lapack::lapack gfortran::gfortran)
    else()
        message(FATAL_ERROR "arpack-ng could not be installed")
    endif()
endif()
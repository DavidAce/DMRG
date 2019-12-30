find_package(arpack-ng HINTS ${CMAKE_INSTALL_PREFIX}/arpack-ng)
if(arpack_ng_LIBRARIES AND arpack_ng_INCLUDE_DIRS)
    add_library(arpack::arpack ${LINK_TYPE} IMPORTED)
    set_target_properties(arpack::arpack PROPERTIES IMPORTED_LOCATION "${arpack_ng_LIBRARIES}")
    target_include_directories(arpack::arpack SYSTEM INTERFACE ${arpack_ng_INCLUDE_DIRS})
    target_link_libraries(arpack::arpack INTERFACE blas::blas lapack::lapack gfortran::gfortran)

endif()

if (NOT TARGET arpack::arpack)
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
            target_link_libraries(arpack::arpack INTERFACE lapacke openblas )
            message(WARNING "Found arpack-ng in system as shared library. "
                            "Make sure to use correct libraries libblas.so.3 and liblapack.so.3 "
                            "by using \n"
                            "sudo update-alternatives --config libblas.so.3-x86_64-linux-gnu \n"
                            "sudo update-alternatives --config liblapack.so.3-x86_64-linux-gnu ")
        endif()
    endif()
endif()


if(NOT TARGET arpack::arpack)
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
    include(ExternalProject)
    include(GNUInstallDirs)
    ExternalProject_Add(external_ARPACK
            GIT_REPOSITORY      https://github.com/opencollab/arpack-ng.git
#            GIT_TAG             master
            GIT_TAG             3.7.0
            GIT_PROGRESS false
            GIT_SHALLOW true
            PREFIX      ${EXTERNAL_BUILD_DIR}/arpack-ng
            INSTALL_DIR ${EXTERNAL_INSTALL_DIR}/arpack-ng
            UPDATE_COMMAND ""
            BUILD_IN_SOURCE 1
            CMAKE_GENERATOR "CodeBlocks - Unix Makefiles"
            CMAKE_ARGS
            -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
            -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
            -DCMAKE_INSTALL_MESSAGE=NEVER #Avoid unnecessary output to console
            -DCMAKE_C_FLAGS=${ARPACK_FLAGS}
            -DCMAKE_Fortran_FLAGS=${ARPACK_FLAGS}
            -DEXAMPLES=ON
            -DCMAKE_BUILD_TYPE=Release
            -DMPI=OFF
            -DINTERFACE64=OFF
            -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
            -DBLAS_LIBRARIES=${EXPANDED_BLAS_GENERATOR}
            -DLAPACK_LIBRARIES=${EXPANDED_LAPACK_GENERATOR}
            DEPENDS blas::blas lapack::lapack gfortran::gfortran
            BUILD_BYPRODUCTS <INSTALL_DIR>/${CMAKE_INSTALL_LIBDIR}/libarpack${ARPACK_SUFFIX}
            )
    ExternalProject_Get_Property(external_ARPACK INSTALL_DIR)
    set(ARPACK_LIBRARIES ${INSTALL_DIR}/${CMAKE_INSTALL_LIBDIR}/libarpack${ARPACK_SUFFIX})
    set(ARPACK_INCLUDE_DIRS ${INSTALL_DIR}/${CMAKE_INSTALL_INCLUDEDIR})

    add_library(arpack::arpack ${LINK_TYPE} IMPORTED)
    set_target_properties(arpack::arpack PROPERTIES IMPORTED_LOCATION "${ARPACK_LIBRARIES}")
    target_include_directories(arpack::arpack SYSTEM INTERFACE ${ARPACK_INCLUDE_DIRS})
    target_link_libraries(arpack::arpack INTERFACE blas::blas lapack::lapack gfortran::gfortran)
    add_dependencies(arpack::arpack external_ARPACK)
endif()
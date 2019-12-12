
unset(ARPACK_LIBRARIES)
find_package(arpack-ng HINTS ${CMAKE_INSTALL_PREFIX}/arpack-ng)
if(arpack_ng_LIBRARIES AND arpack_ng_INCLUDE_DIRS)
    add_library(arpack INTERFACE IMPORTED)
    target_link_libraries(arpack INTERFACE ${arpack_ng_LIBRARIES} blas lapack gfortran)
    target_include_directories(arpack SYSTEM INTERFACE ${arpack_ng_INCLUDE_DIRS})
endif()

if (NOT TARGET arpack)
    message(STATUS "Searching for Arpack-ng")
    find_library(ARPACK_LIBRARIES
            NAMES libarpack${CUSTOM_SUFFIX} arpack
            PATH_SUFFIXES lib lib32 lib64 x86_64-linux-gnu lib/x86_64-linux-gnu
            HINTS
                $ENV{ARPACK_DIR}
                ${CMAKE_INSTALL_PREFIX}/arpack-ng
                $ENV{CONDA_PREFIX}
            PATHS
                $ENV{EBROOTARPACKMINNG}
            PATH_SUFFIXES arpack/lib arpack-ng/lib lib lib32 lib64
            )
    if(NOT ARPACK_LIBRARIES)
        message(STATUS "Searching for Arpack-ng - failed")
    else()
        message(STATUS "Searching for Arpack-ng - Success: ${ARPACK_LIBRARIES}")
        add_library(arpack INTERFACE IMPORTED)
        target_link_libraries(arpack INTERFACE ${ARPACK_LIBRARIES} blas lapack gfortran)
    endif()
endif()


if(NOT TARGET arpack)
    message(STATUS "Arpack-ng will be installed into ${CMAKE_INSTALL_PREFIX}/arpack-ng")

    #####################################################################
    ### Prepare lists with generator expressions, replacing all semicolons.
    ### Otherwise, passing raw lists results  in only the first element
    ### of the list to be passed.
    ####################################################################
    include(${PROJECT_SOURCE_DIR}/cmake-modules/getExpandedTarget.cmake)
    include(${PROJECT_SOURCE_DIR}/cmake-modules/filterTarget.cmake)
    add_library(arpack-aux-target INTERFACE IMPORTED)
    target_link_libraries(arpack-aux-target INTERFACE blas lapack)
    remove_pthread(arpack-aux-target)
    remove_duplicates(arpack-aux-target)
    expand_target_libs(arpack-aux-target AUX_LIBRARIES)
    message("AUX_LIBRARIES   : ${AUX_LIBRARIES}")
    string (REPLACE ";" "$<SEMICOLON>" AUX_LIBRARIES_GENERATOR     "${AUX_LIBRARIES}")

    message("BLAS_LIBRARIES_GENERATOR   : ${AUX_LIBRARIES_GENERATOR}")
    ####################################################################
    set(ARPACK_FLAGS "-w -m64 -fPIC")
    include(ExternalProject)
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
            -DBLAS_LIBRARIES=${AUX_LIBRARIES_GENERATOR}
            -DLAPACK_LIBRARIES=${AUX_LIBRARIES_GENERATOR}
            DEPENDS blas lapack gfortran
            )
    ExternalProject_Get_Property(external_ARPACK INSTALL_DIR)
    include(GNUInstallDirs)
    set(ARPACK_LIBRARIES ${INSTALL_DIR}/${CMAKE_INSTALL_LIBDIR}/libarpack${CUSTOM_SUFFIX})
    set(ARPACK_INCLUDE_DIRS ${INSTALL_DIR}/${CMAKE_INSTALL_INCLUDEDIR})
    add_library(arpack INTERFACE IMPORTED)
    add_dependencies(arpack external_ARPACK)
    target_link_libraries(arpack INTERFACE  ${ARPACK_LIBRARIES} blas lapack gfortran)
    target_include_directories(arpack SYSTEM INTERFACE  ${ARPACK_INCLUDE_DIRS})
endif()
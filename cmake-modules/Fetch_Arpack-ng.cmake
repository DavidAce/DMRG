
unset(ARPACK_LIBRARIES)

if(BUILD_SHARED_LIBS)
    set(ARPACK_LIBRARY_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
else()
    set(ARPACK_LIBRARY_SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()

if (NOT ARPACK_LIBRARIES)
message(STATUS "Searching for Arpack-ng in system")
find_library(ARPACK_LIBRARIES
        NAMES libarpack${ARPACK_LIBRARY_SUFFIX}
        PATH_SUFFIXES lib lib32 lib64 x86_64-linux-gnu lib/x86_64-linux-gnu
        PATHS /usr /usr/local
        NO_DEFAULT_PATH
    )
    if(NOT ARPACK_LIBRARIES)
    message(STATUS "Searching for Arpack-ng in system - failed")
    else()
    message(STATUS "Searching for Arpack-ng in system - Success: ${ARPACK_LIBRARIES}")
    endif()
endif()

if (NOT ARPACK_LIBRARIES)
    # Try finding arpack as module library
    message(STATUS "Searching for Arpack-ng in module")
    find_library(ARPACK_LIBRARIES
            NAMES libarpack${ARPACK_LIBRARY_SUFFIX} arpack
            PATH_SUFFIXES lib lib32 lib64
            PATHS
            $ENV{EBROOTARPACKMINNG}
            $ENV{ARPACK_DIR}
            NO_DEFAULT_PATH
            )
    if(NOT ARPACK_LIBRARIES)
        message(STATUS "Searching for Arpack-ng in module - failed")
    else()
        message(STATUS "Searching for Arpack-ng in module - Success: ${ARPACK_LIBRARIES}")
    endif()
endif()



if(ARPACK_LIBRARIES)
    add_library(arpack INTERFACE)
    target_link_libraries(arpack INTERFACE ${ARPACK_LIBRARIES} blas lapack gfortran)
    add_dependencies(arpack blas lapack gfortran)
else()
    message(STATUS "Arpack-ng will be installed into ${EXTERNAL_INSTALL_DIR}/arpack-ng on first build.")

    #####################################################################
    ### Prepare lists with generator expressions, replacing all semicolons.
    ### Otherwise, passing raw lists results  in only the first element
    ### of the list to be passed.
    ####################################################################
#    get_target_property(BLAS_LIBRARIES      blas    INTERFACE_LINK_LIBRARIES)
#    get_target_property(LAPACK_LIBRARIES    lapack  INTERFACE_LINK_LIBRARIES)
    include(${PROJECT_SOURCE_DIR}/cmake-modules/getExpandedTarget.cmake)
    expand_target_libs(blas BLAS_LIBRARIES)
    expand_target_libs(lapack LAPACK_LIBRARIES)

    string (REPLACE ";" "$<SEMICOLON>" BLAS_LIBRARIES_GENERATOR     "${BLAS_LIBRARIES}")
    string (REPLACE ";" "$<SEMICOLON>" LAPACK_LIBRARIES_GENERATOR   "${LAPACK_LIBRARIES}")
#    string (REPLACE ";" "$<SEMICOLON>" FC_LDLAGS_GENERATOR          "${FC_LDLAGS}")
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
            -DCMAKE_CXX_STANDARD=17
            -DCMAKE_CXX_STANDARD_REQUIRED:BOOL=ON
            -DCMAKE_CXX_EXTENSIONS:BOOL=OFF
            -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
            -DCMAKE_INSTALL_MESSAGE=NEVER #Avoid unnecessary output to console
            -DCMAKE_C_FLAGS=${ARPACK_FLAGS}
            -DCMAKE_Fortran_FLAGS=${ARPACK_FLAGS}
            -DEXAMPLES=ON
            -DCMAKE_BUILD_TYPE=Release
            -DMPI=OFF
            -DINTERFACE64=OFF
            -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
            -DBLAS_LIBRARIES=${BLAS_LIBRARIES_GENERATOR}
            -DLAPACK_LIBRARIES=${LAPACK_LIBRARIES_GENERATOR}
#            -DEXTRA_LDLAGS=${FC_LDLAGS_GENERATOR}
            DEPENDS blas lapack gfortran
            )
    ExternalProject_Get_Property(external_ARPACK INSTALL_DIR)
    include(GNUInstallDirs)
    set(ARPACK_LIBRARIES ${INSTALL_DIR}/${CMAKE_INSTALL_LIBDIR}/libarpack${ARPACK_LIBRARY_SUFFIX})
    set(ARPACK_INCLUDE_DIRS ${INSTALL_DIR}/${CMAKE_INSTALL_INCLUDEDIR})
    add_library(arpack INTERFACE)
    add_dependencies(arpack external_ARPACK)
    add_dependencies(arpack blas lapack gfortran)
#    message(STATUS "ARPACK_LIBRARIES   :  ${ARPACK_LIBRARIES}")
#    message(STATUS "ARPACK_INCLUDE_DIRS:  ${ARPACK_INCLUDE_DIRS}")
    target_link_libraries(arpack INTERFACE  ${ARPACK_LIBRARIES} blas lapack )
    target_include_directories(arpack SYSTEM INTERFACE  ${ARPACK_INCLUDE_DIRS})
endif()
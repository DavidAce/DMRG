
unset(ARPACK_LIBRARIES)

if(BUILD_SHARED_LIBS)
    set(ARPACK_LIBRARY_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
else()
    set(ARPACK_LIBRARY_SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()

message(STATUS "SEARCHING FOR ARPACK IN SYSTEM")
find_library(ARPACK_LIBRARIES
        NAMES libarpack${ARPACK_LIBRARY_SUFFIX}
        PATHS
            /usr/lib/x86_64-linux-gnu/
            /usr/lib
        NO_DEFAULT_PATH
    )

if (NOT ARPACK_LIBRARIES)
    # Try finding arpack as module library
    message(STATUS "ARPACK NOT FOUND IN SYSTEM. SEARCH RETURNED:  ${ARPACK_LIBRARIES}")
    message(STATUS "SEARCHING FOR ARPACK IN LOADED MODULES")
    find_library(ARPACK_LIBRARIES
            NAMES libarpack${ARPACK_LIBRARY_SUFFIX} arpack
            PATHS "$ENV{ARPACK_DIR}/lib" "$ENV{ARPACK_DIR}/lib64"
            NO_DEFAULT_PATH
            )
    find_path(ARPACK_INCLUDE_DIRS
            NAMES arpack
            PATHS "$ENV{ARPACK_DIR}/include"
            NO_DEFAULT_PATH
            )
else()
    set(ARPACK_INCLUDE_DIRS TRUE)
endif()



if(ARPACK_LIBRARIES AND ARPACK_INCLUDE_DIRS)
    message(STATUS "ARPACK found in system:   ${ARPACK_LIBRARIES}")
    add_library(arpack INTERFACE)
    target_link_libraries(arpack INTERFACE ${ARPACK_LIBRARIES} blas lapack gfortran)
    target_include_directories(arpack INTERFACE ${ARPACK_INCLUDE_DIRS})
    add_dependencies(arpack blas lapack gfortran)
else()
    message(STATUS "Arpack-ng will be installed into ${INSTALL_DIRECTORY}/arpack-ng on first build.")

    #####################################################################
    ### Prepare lists with generator expressions, replacing all semicolons.
    ### Otherwise, passing raw lists results  in only the first element
    ### of the list to be passed.
    ####################################################################

    string (REPLACE ";" "$<SEMICOLON>" BLAS_LIBRARIES_GENERATOR     "${BLAS_LIBRARIES}")
    string (REPLACE ";" "$<SEMICOLON>" LAPACK_LIBRARIES_GENERATOR   "${LAPACK_LIBRARIES}")
    string (REPLACE ";" "$<SEMICOLON>" FC_LDLAGS_GENERATOR          "${FC_LDLAGS}")
    ####################################################################

    include(ExternalProject)
    ExternalProject_Add(external_ARPACK
            GIT_REPOSITORY      https://github.com/opencollab/arpack-ng.git
#            GIT_TAG             master
            GIT_TAG             3.7.0
            PREFIX      ${BUILD_DIRECTORY}/arpack-ng
            INSTALL_DIR ${INSTALL_DIRECTORY}/arpack-ng
            UPDATE_COMMAND ""
            BUILD_IN_SOURCE 1
            CMAKE_ARGS
            -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
            -DCMAKE_INSTALL_MESSAGE=NEVER #Avoid unnecessary output to console
            -DCMAKE_C_FLAGS=-w -m64
            -DCMAKE_Fortran_FLAGS=-w -m64
            -DEXAMPLES=ON
            -DCMAKE_BUILD_TYPE=Release
            -DMPI=OFF
            -DINTERFACE64=OFF
            -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
            -DBLAS_LIBRARIES=${BLAS_LIBRARIES_GENERATOR}
            -DLAPACK_LIBRARIES=${LAPACK_LIBRARIES_GENERATOR}
            -DEXTRA_LDLAGS=${FC_LDLAGS_GENERATOR}
            DEPENDS blas lapack gfortran
            )
    ExternalProject_Get_Property(external_ARPACK INSTALL_DIR)
    include(GNUInstallDirs)
    set(ARPACK_LIBRARIES ${INSTALL_DIR}/${CMAKE_INSTALL_LIBDIR}/libarpack${ARPACK_LIBRARY_SUFFIX})
    set(ARPACK_INCLUDE_DIRS ${INSTALL_DIR}/${CMAKE_INSTALL_INCLUDEDIR})
    add_library(arpack INTERFACE)
    add_dependencies(arpack external_ARPACK)
    add_dependencies(arpack blas lapack gfortran)
    message(STATUS "ARPACK_LIBRARIES   :  ${ARPACK_LIBRARIES}")
    message(STATUS "ARPACK_INCLUDE_DIRS:  ${ARPACK_INCLUDE_DIRS}")
    target_link_libraries(arpack INTERFACE  ${ARPACK_LIBRARIES} blas lapack )
    target_include_directories(arpack INTERFACE ${ARPACK_INCLUDE_DIRS})
endif()



message(STATUS "SEARCHING FOR ARPACK IN SYSTEM...")
find_library(ARPACK_LIBRARIES
        NAMES libarpack${CUSTOM_SUFFIX}
        PATHS /usr/lib/x86_64-linux-gnu/
        NO_DEFAULT_PATH
        )

if (NOT ARPACK_LIBRARIES)
    # Try finding arpack as module library
    message(STATUS "SEARCHING FOR ARPACK IN LOADED MODULES")

    find_library(ARPACK_LIBRARIES
            NAMES libarpack${CUSTOM_SUFFIX} arpack
            PATHS "$ENV{ARPACK_DIR}/lib" "$ENV{ARPACK_DIR}/lib64"
            NO_DEFAULT_PATH
            )
    find_path(ARPACK_INCLUDE_DIRS
            NAMES arpack
            PATHS "$ENV{ARPACK_DIR}/include"
            NO_DEFAULT_PATH
            )
endif()



if(ARPACK_LIBRARIES AND ARPACK_INCLUDE_DIRS)
    message(STATUS "ARPACK found in system:   ${ARPACK_LIBRARIES}")
    add_library(arpack INTERFACE)
    set_target_properties(arpack
            PROPERTIES
            INTERFACE_LINK_LIBRARIES      "${ARPACK_LIBRARIES};blas;lapack;gfortran"
            INTERFACE_INCLUDE_DIRECTORIES "${ARPACK_INCLUDE_DIRS}"
            #            INTERFACE_LINK_OPTIONS        "${PTHREAD_LIBRARY}"
            )

    add_dependencies(arpack blas lapack gfortran)
    return()
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
            GIT_REPOSITORY      git@github.com:opencollab/arpack-ng.git
#            GIT_TAG             master
            GIT_TAG             3.6.3
            PREFIX      ${BUILD_DIRECTORY}/arpack-ng
            INSTALL_DIR ${INSTALL_DIRECTORY}/arpack-ng
            UPDATE_COMMAND ""
            BUILD_IN_SOURCE 1
#            CONFIGURE_COMMAND
#                ./bootstrap && export INTERFACE64=1 &&
#                ./configure INTERFACE64=1
#                            --prefix=<INSTALL_DIR>
#                            --enable-silent-rules
#                            --with-blas=${BLAS_LIBRARIES_GENERATOR}
#                            --with-lapack=${LAPACK_LIBRARIES_GENERATOR}
#
#            BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} && ${CMAKE_MAKE_PROGRAM} check
#            INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install

            CMAKE_ARGS
            -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
            -DCMAKE_INSTALL_MESSAGE=NEVER #Avoid unnecessary output to console
            -DCMAKE_C_FLAGS=-w -m64 -fPIC
            -DCMAKE_Fortran_FLAGS=-w -m64 -fPIC
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
    find_library(ARPACK_LIBRARIES
            NAMES libarpack${CUSTOM_SUFFIX}
            PATHS  ${INSTALL_DIR}/lib  ${INSTALL_DIR}/lib64
            )
    if(NOT ARPACK_LIBRARIES)
        set(ARPACK_LIBRARIES ${INSTALL_DIR}/lib/libarpack${CUSTOM_SUFFIX})
    endif()
    set(ARPACK_INCLUDE_DIRS ${INSTALL_DIR}/include)
    message("THIS IS THE PATH MAKING ALL THE FUSS: ${ARPACK_INCLUDE_DIRS}")
    add_library(arpack INTERFACE)
    set_target_properties(arpack
            PROPERTIES
            INTERFACE_LINK_LIBRARIES      "${ARPACK_LIBRARIES};blas;lapack;gfortran"
            INTERFACE_INCLUDE_DIRECTORIES "${ARPACK_INCLUDE_DIRS}"
            )

    add_dependencies(arpack external_ARPACK blas lapack gfortran)
    #target_link_libraries(${PROJECT_NAME} PRIVATE arpack)
endif()
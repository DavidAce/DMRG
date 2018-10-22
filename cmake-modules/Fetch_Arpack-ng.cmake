


message(STATUS "SEARCHING FOR ARPACK IN SYSTEM...")
find_library(ARPACK_LIBRARIES
        NAMES libarpack${CUSTOM_SUFFIX}
        PATHS /usr/lib/x86_64-linux-gnu/
        )

if (NOT ARPACK_LIBRARIES)
    # Try finding arpack as module library
    message(STATUS "SEARCHING FOR ARPACK IN LOADED MODULES")

    find_library(ARPACK_LIBRARIES
            NAMES arpack
            PATHS "$ENV{ARPACK_DIR}/lib"
            )
    find_path(ARPACK_INCLUDE_DIRS
            NAMES arpack
            PATHS "$ENV{ARPACK_DIR}/include"
            )
endif()



if(ARPACK_LIBRARIES)
    message(STATUS "ARPACK found in system:   ${ARPACK_LIBRARIES}")
    add_library(arpack UNKNOWN IMPORTED)
    set_target_properties(arpack
            PROPERTIES
            IMPORTED_LOCATION "${ARPACK_LIBRARIES}"
            INTERFACE_LINK_LIBRARIES "blas;lapack;gfortran;-lpthread"
            INTERFACE_INCLUDE_DIRECTORIES "${ARPACK_INCLUDE_DIRS}"
            INTERFACE_LINK_FLAGS          "${OpenMP_CXX_FLAGS}"
            )
    target_link_libraries(${PROJECT_NAME} PRIVATE arpack)
    target_include_directories(${PROJECT_NAME} PRIVATE ${ARPACK_INCLUDE_DIRS})
#    target_link_libraries(${PROJECT_NAME} PRIVATE -lpthread)

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
    if(${STATIC_BUILD})
        set(ARPACK_SHARED OFF)
    else()
        set(ARPACK_SHARED ON)
    endif()

    include(ExternalProject)
    ExternalProject_Add(library_ARPACK
            GIT_REPOSITORY      https://github.com/opencollab/arpack-ng.git
            GIT_TAG             master # You need to do shared library linking with blas and lapack for this to work, otherwise the examples will fail due to missing -lpthread
#            GIT_TAG             3.5.0 # Latest version has problems with fortran linking. so stick with this version instead.
            PREFIX              "${INSTALL_DIRECTORY}/arpack-ng"
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
            -DBUILD_SHARED_LIBS=${ARPACK_SHARED}
            -DBLAS_LIBRARIES=${BLAS_LIBRARIES_GENERATOR};
            -DLAPACK_LIBRARIES=${LAPACK_LIBRARIES_GENERATOR}
            -DEXTRA_LDLAGS=${FC_LDLAGS_GENERATOR}
            DEPENDS blas lapack gfortran
            )
    ExternalProject_Get_Property(library_ARPACK INSTALL_DIR)
    set(ARPACK_INCLUDE_DIRS ${INSTALL_DIR}/include)
    add_library(arpack STATIC IMPORTED)
    set_target_properties(arpack
            PROPERTIES
            IMPORTED_LOCATION "${INSTALL_DIR}/lib/libarpack${CUSTOM_SUFFIX}"
            INTERFACE_LINK_LIBRARIES "blas;lapack;gfortran"
            INTERFACE_LINK_FLAGS            "-lpthread"
            )
#            INTERFACE_LINK_FLAGS      )

    add_dependencies(arpack library_ARPACK blas lapack gfortran )
    target_link_libraries(${PROJECT_NAME} PRIVATE arpack -lpthread -m64)
    target_include_directories(${PROJECT_NAME} PRIVATE ${ARPACK_INCLUDE_DIRS})
    #For convenience, define these variables
    get_target_property(ARPACK_LIBRARIES arpack IMPORTED_LOCATION)
endif()
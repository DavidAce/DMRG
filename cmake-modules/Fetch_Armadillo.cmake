
# Try finding Armadillo in system
find_package(Armadillo)

if (NOT ARMADILLO_FOUND)
# Try finding armadillo as module library
find_library(ARMADILLO_LIBRARIES
        NAMES armadillo
        PATHS "$ENV{ARMADILLO_DIR}/lib"
        )
find_path(ARMADILLO_INCLUDE_DIRS
        NAMES armadillo
        PATHS "$ENV{ARMADILLO_DIR}/include"
        )
endif()
if(ARMADILLO_FOUND OR ARMADILLO_LIBRARIES)
    message(STATUS "ARMADILLO found in system: ${ARMADILLO_LIBRARIES}")
    add_library(armadillo STATIC IMPORTED)
#else()
#    # Else try finding a previously installed Armadillo in local folder
#    find_library(ARMADILLO_LIBRARIES NAMES armadillo libarmadillo.so PATHS "${INSTALL_DIRECTORY}/armadillo/lib" )
#    if(ARMADILLO_LIBRARIES)
#        message(STATUS "ARMADILLO already installed: ${ARMADILLO_LIBRARIES}")
#        set(ARMADILLO_INCLUDE_DIRS ${INSTALL_DIRECTORY}/armadillo/include)
#        add_library(armadillo STATIC IMPORTED)
#    endif()
endif()

# Install from source if not found
if(NOT ARMADILLO_LIBRARIES)
    message(STATUS "ARMADILLO will be installed into ${INSTALL_DIRECTORY}/armadillo on first build.")
    #####################################################################
    ### Prepare lists with generator expressions, replacing all semicolons.
    ### Otherwise, passing raw lists results  in only the first element
    ### of the list to be passed.
    ####################################################################
    string (REPLACE ";" "$<SEMICOLON>" BLAS_LIBRARIES_GENERATOR     "${BLAS_LIBRARIES};${EXTRA_LDLAGS}")
    string (REPLACE ";" "$<SEMICOLON>" LAPACK_LIBRARIES_GENERATOR   "${LAPACK_LIBRARIES}")
    string (REPLACE ";" "$<SEMICOLON>" EXTRA_LDLAGS_GENERATOR       "${EXTRA_LDLAGS}")
    ####################################################################

    include(ExternalProject)
    ExternalProject_Add(library_ARMADILLO
            GIT_REPOSITORY      https://gitlab.com/conradsnicta/armadillo-code.git
            GIT_TAG             9.100.x
            PREFIX              "${INSTALL_DIRECTORY}/armadillo"
#            UPDATE_DISCONNECTED 0
            UPDATE_COMMAND ""
            TEST_COMMAND ""
            CMAKE_ARGS
            -DBUILD_SHARED_LIBS=OFF
            -DCMAKE_SYSTEM_LIBRARY_PATH=${INSTALL_DIRECTORY}/OpenBLAS/lib
            -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
            -DARMA_USE_LAPACK=true
            -DARMA_USE_BLAS=true
            -DBLAS_LIBRARIES=${BLAS_LIBRARIES_GENERATOR}
            -DLAPACK_LIBRARY=${LAPACK_LIBRARIES_GENERATOR}
            -DLAPACK_NAMES=${LAPACK_LIBRARIES_GENERATOR}
            -DDETECT_HDF5=OFF
            DEPENDS arpack blas lapack gfortran
            )

    ExternalProject_Get_Property(library_ARMADILLO INSTALL_DIR)
    add_library(armadillo           STATIC IMPORTED)
    add_dependencies(armadillo      library_ARMADILLO arpack blas lapack)
    set(ARMADILLO_LIBRARIES         ${INSTALL_DIR}/lib/libarmadillo${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(ARMADILLO_INCLUDE_DIRS      ${INSTALL_DIR}/include)


endif()

set_target_properties(armadillo PROPERTIES
        IMPORTED_LOCATION        "${ARMADILLO_LIBRARIES}"
        INTERFACE_LINK_LIBRARIES "blas;lapack;gfortran"
        INCLUDE_DIRECTORIES      "${ARMADILLO_INCLUDE_DIRS}")


target_link_libraries(      ${PROJECT_NAME} PRIVATE armadillo)
target_include_directories( ${PROJECT_NAME} PRIVATE ${ARMADILLO_INCLUDE_DIRS})
target_compile_definitions( ${PROJECT_NAME} PRIVATE -DARMA_NO_DEBUG)

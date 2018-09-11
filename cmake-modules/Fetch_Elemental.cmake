
enable_language(Fortran)
include(cmake-modules/FindGFortran.cmake)
find_file(ELEMENTAL_FOUND PATHS ${INSTALL_DIRECTORY}/elemental/lib NAMES libEl.so NO_DEFAULT_PATH)
if(ELEMENTAL_FOUND)
    add_library(elemental           SHARED IMPORTED)
else()
    message(STATUS "ELEMENTAL will be installed into ${INSTALL_DIRECTORY}/elemental on first build.")

    #####################################################################
    ### Prepare lists with generator expressions, replacing all semicolons.
    ### Otherwise, passing raw lists results  in only the first element
    ### of the list to be passed.
    ####################################################################
    string (REPLACE ";" "$<SEMICOLON>" BLAS_LIBRARIES_GENERATOR     "${BLAS_LIBRARIES}")
    string (REPLACE ";" "$<SEMICOLON>" LAPACK_LIBRARIES_GENERATOR   "${LAPACK_LIBRARIES}")
    string (REPLACE ";" "$<SEMICOLON>" EXTRA_LDLAGS_GENERATOR       "${EXTRA_LDLAGS}")
#    set(BLA_VENDOR "OpenBLAS")
    ####################################################################

    include(ExternalProject)
    ExternalProject_Add(library_ELEMENTAL
            GIT_REPOSITORY      https://github.com/elemental/Elemental.git
            GIT_TAG             v0.87.7
            PREFIX              "${INSTALL_DIRECTORY}/elemental"
            UPDATE_DISCONNECTED 1
            TEST_COMMAND ""

            CMAKE_ARGS
                -DEL_C_INTERFACE=OFF
                -DEL_HYBRID=OFF
                -DEL_DISABLE_PARMETIS=ON
                -DGFORTRAN_LIB=${GFORTRAN_LIB}
                -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
                -DCMAKE_BUILD_TYPE=Release
                -DBLAS_LIBRARIES=${BLAS_LIBRARIES_GENERATOR}
                -DLAPACK_LIBRARIES=${LAPACK_LIBRARIES_GENERATOR}
                -DBLAS_PATH=${BLAS_LIBRARIES_GENERATOR}
                -DLAPACK_PATH=${LAPACK_LIBRARIES_GENERATOR}
                -DMATH_LIBS=${BLAS_LIBRARIES_GENERATOR}
                -DMATH_PATH=${BLAS_LIBRARIES_GENERATOR}
                -DCMAKE_REQUIRED_LIBRARIES=${BLAS_LIBRARIES_GENERATOR}

            DEPENDS blas lapack
            )

    add_library(elemental           SHARED IMPORTED)
    add_dependencies(elemental      library_ELEMENTAL blas lapack)

endif()
set(INSTALL_DIR ${INSTALL_DIRECTORY}/elemental)
set(ELEMENTAL_LIBRARIES                 ${INSTALL_DIR}/lib/libEl${CMAKE_SHARED_LIBRARY_SUFFIX})
set(ELEMENTAL_LINK_LIBRARIES            ${INSTALL_DIR}/lib/libElSuiteSparse${CMAKE_SHARED_LIBRARY_SUFFIX})
list(APPEND ELEMENTAL_LINK_LIBRARIES    ${INSTALL_DIR}/lib/libmetis${CMAKE_SHARED_LIBRARY_SUFFIX})
#    list(APPEND ELEMENTAL_LINK_LIBRARIES    ${INSTALL_DIR}/lib/libparmetis${CMAKE_SHARED_LIBRARY_SUFFIX})
list(APPEND ELEMENTAL_LINK_LIBRARIES    ${INSTALL_DIR}/lib/libpmrrr${CMAKE_SHARED_LIBRARY_SUFFIX})
set(ELEMENTAL_INCLUDE_DIRS      ${INSTALL_DIR}/include)


###################
####  Link MPI  ###
###################
find_package(MPI REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE ${MPI_LIBRARIES})
set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
target_include_directories(${PROJECT_NAME} PRIVATE ${MPI_INCLUDE_PATH})
#target_link_libraries(${PROJECT_NAME} PRIVATE ${MPI_LIBRARIES})
#set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
#set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
#target_include_directories(${PROJECT_NAME} PRIVATE ${MPI_INCLUDE_PATH})
#


set_target_properties(elemental PROPERTIES
        IMPORTED_LOCATION "${ELEMENTAL_LIBRARIES}"
        INTERFACE_LINK_LIBRARIES "${ELEMENTAL_LINK_LIBRARIES}"
        INCLUDE_DIRECTORIES "${ELEMENTAL_INCLUDE_DIRS}")


target_link_libraries(${PROJECT_NAME} PRIVATE elemental)

# Add SYSTEM flag to suppress warnings
target_include_directories(${PROJECT_NAME} SYSTEM PRIVATE ${ELEMENTAL_INCLUDE_DIRS})


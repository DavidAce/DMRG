
# Try finding in system. While avoiding Python anaconda's version.
# Usually there is an executable present in system, "h5c++"

find_package(Armadillo)
if(ARMADILLO_FOUND)
    add_library(armadillo SHARED IMPORTED)
endif()


if(NOT ARMADILLO_FOUND)
    message(STATUS "ARMADILLO will be installed into ${INSTALL_DIRECTORY}/armadillo on first build.")

    include(ExternalProject)
    ExternalProject_Add(library_ARMADILLO
            GIT_REPOSITORY      https://gitlab.com/conradsnicta/armadillo-code.git
            GIT_TAG             8.600.x
            PREFIX              "${INSTALL_DIRECTORY}/armadillo"
            UPDATE_DISCONNECTED 1
            TEST_COMMAND ""
            CMAKE_ARGS
            -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
            -DDETECT_HDF5=OFF
            DEPENDS arpack blas lapack
            )

    ExternalProject_Get_Property(library_ARMADILLO INSTALL_DIR)
    add_library(armadillo           SHARED IMPORTED)
    add_dependencies(armadillo     library_ARMADILLO arpack blas lapack)
    set(ARMADILLO_LIBRARIES        ${INSTALL_DIR}/lib/libarmadillo${CMAKE_SHARED_LIBRARY_SUFFIX})
    set(ARMADILLO_INCLUDE_DIRS      ${INSTALL_DIR}/include)


endif()

set_target_properties(armadillo PROPERTIES
        IMPORTED_LOCATION "${ARMADILLO_LIBRARIES}"
        INTERFACE_LINK_LIBRARIES "${BLAS_LIBRARIES};${LAPACK_LIBRARIES}"
        INCLUDE_DIRECTORIES "${ARMADILLO_INCLUDE_DIRS}")


target_link_libraries(${PROJECT_NAME} PRIVATE armadillo)
target_include_directories(${PROJECT_NAME} PRIVATE ${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PRIVATE ${BLAS_LIBRARIES})
target_compile_definitions(${PROJECT_NAME} PRIVATE -DARMA_NO_DEBUG)

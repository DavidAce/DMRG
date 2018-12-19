
# Try finding Armadillo in system
find_package(Armadillo)

if (NOT ARMADILLO_FOUND)
# Try finding armadillo as module library
find_path(ARMADILLO_INCLUDE_DIRS
        NAMES armadillo
        PATHS "$ENV{ARMADILLO_DIR}/include"
        )
endif()
if(ARMADILLO_INCLUDE_DIRS)
    message(STATUS "ARMADILLO found in system: ${ARMADILLO_INCLUDE_DIRS}")
    add_library(armadillo INTERFACE)
endif()

# Install from source if not found
if(NOT ARMADILLO_INCLUDE_DIRS)
    message(STATUS "ARMADILLO will be installed into ${INSTALL_DIRECTORY}/armadillo on first build.")
    include(ExternalProject)
    ExternalProject_Add(library_ARMADILLO
            GIT_REPOSITORY      https://gitlab.com/conradsnicta/armadillo-code.git
            GIT_TAG             9.200.x
            PREFIX              "${INSTALL_DIRECTORY}/armadillo"
#            UPDATE_DISCONNECTED 0
            UPDATE_COMMAND ""
            CONFIGURE_COMMAND sed -i "s/^set(ARMA_USE_WRAPPER true)/set(ARMA_USE_WRAPPER false)/" <SOURCE_DIR>/CMakeLists.txt
            COMMAND cmake
                -DBUILD_SHARED_LIBS=OFF
                -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
                -DDETECT_HDF5=OFF
                <SOURCE_DIR>
            TEST_COMMAND ""

            DEPENDS arpack blas lapack gfortran
            )

    ExternalProject_Get_Property(library_ARMADILLO INSTALL_DIR)
    add_library(armadillo           INTERFACE)
    add_dependencies(armadillo      library_ARMADILLO arpack blas lapack eigen3)
#    set(ARMADILLO_LIBRARIES         ${INSTALL_DIR}/lib/libarmadillo${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(ARMADILLO_INCLUDE_DIRS      ${INSTALL_DIR}/include)

endif()

set_target_properties(armadillo PROPERTIES
        INTERFACE_LINK_LIBRARIES         "blas;lapack;gfortran;eigen3"
        INTERFACE_COMPILE_OPTIONS        "-DARMA_DONT_USE_WRAPPER;-DARMA_NO_DEBUG"
        INTERFACE_INCLUDE_DIRECTORIES    "${ARMADILLO_INCLUDE_DIRS}")


target_link_libraries(      ${PROJECT_NAME} PRIVATE armadillo)

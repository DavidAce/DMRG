
find_package(h5pp NO_DEFAULT_PATH HINTS ${INSTALL_DIRECTORY}/h5pp $ENV{H5PP_DIR})

if(h5pp_FOUND)
    message(STATUS "h5pp FOUND IN SYSTEM: ${H5PP_INCLUDE_DIR}")
#    add_library(h5pp ALIAS h5pp::h5pp)
else()
    message(STATUS "h5pp will be installed into ${INSTALL_DIRECTORY}/h5pp on first build.")

    include(ExternalProject)
    ExternalProject_Add(external_H5PP
            GIT_REPOSITORY git@github.com:DavidAce/h5pp.git
            GIT_TAG master
            GIT_PROGRESS 1
            PREFIX "${INSTALL_DIRECTORY}/h5pp"
            UPDATE_COMMAND ""
#            CMAKE_ARGS
#            -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
#            -DBUILD_SHARED_LIBS:BOOL=OFF
#            -DDOWNLOAD_ALL:BOOL=OFF
#            -DENABLE_TESTS:BOOL=ON
#            -DBUILD_EXAMPLES:BOOL=OFF
#            -DEigen3_DIR:PATH=${Eigen3_DIR}
#            -DHDF5_DIR:PATH=${HDF5_DIR}
#            -DHDF5_ROOT:PATH=${HDF5_DIR}
#            -Dspdlog_DIR:PATH=${spdlog_DIR}
            CMAKE_ARGS
            -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
            -DMARCH=${MARCH}
            -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
            -DHDF5_DIR:PATH=${HDF5_DIR}
            -DHDF5_ROOT:PATH=${HDF5_DIR}
            -DEigen3_DIR:PATH=${Eigen3_DIR}
            -DEIGEN3_ROOT_DIR:PATH=${EIGEN3_ROOT_DIR}
            -DEIGEN3_INCLUDE_DIR:PATH=${EIGEN3_INCLUDE_DIR}
            -Dspdlog_DIR:PATH=${spdlog_DIR}
            -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
            DEPENDS hdf5 spdlog::spdlog Eigen3::Eigen
            )

    ExternalProject_Get_Property(external_H5PP INSTALL_DIR)
    add_library(h5pp INTERFACE)
    add_library(h5pp::h5pp ALIAS h5pp)

    set(H5PP_INCLUDE_DIR ${INSTALL_DIR}/include/h5pp)
    add_dependencies(h5pp external_H5PP)
    target_include_directories(h5pp INTERFACE ${H5PP_INCLUDE_DIR})
    target_link_libraries(h5pp INTERFACE Eigen3 hdf5 spdlog)
endif()

#This script is supposed to be run from an execute_script() command such as:
#

#        find_package(HDF5 COMPONENTS CXX)
#        if (HDF5_FOUND AND NOT "${HDF5_LIBRARIES}" MATCHES "anaconda")
#            add_definitions    (${HDF5_DEFINITIONS})
#            include_directories(${HDF5_INCLUDE_DIRS})
#            message("FOUND NATIVE HDF5:   ${HDF5_LIBRARIES}")
#        else()
#            message("DOWNLOADING HDF5...")
#            execute_process(
#                    COMMAND ${CMAKE_COMMAND} -E make_directory build
#                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/hdf5)
#            execute_process(
#                    COMMAND ${CMAKE_COMMAND}
#                    -DINSTALL_DIRECTORY:PATH=${PROJECT_SOURCE_DIR}/libs
#                    -G ${CMAKE_GENERATOR} ../
#                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/hdf5/build )
#            execute_process(
#                    COMMAND ${CMAKE_COMMAND} --build .
#                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/hdf5/build )
#            include(${PROJECT_SOURCE_DIR}/libs/hdf5/FindHDF5.cmake)
#            include_directories(${HDF5_INC_DIR})
#            message("SUCCESSFULLY INSTALLED HDF5:   ${HDF5_LIBRARIES}")
#        endif()
# It will generate a file called FindHDF5.cmake in the install directory ${PROJECT_SOURCE_DIR}/libs/hdf5/FindHDF5.cmake
#



cmake_minimum_required(VERSION 3.7)
project(HDF5-download NONE)
include(ExternalProject)
ExternalProject_Add(HDF5
        URL      https://fossies.org/linux/misc/hdf5-1.10.1.tar.bz2
        PREFIX              "${INSTALL_DIRECTORY}/hdf5"
        CMAKE_ARGS
        TEST_COMMAND ""
        CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
        -DBUILD_SHARED_LIBS:BOOL=OFF
        -DBUILD_TESTING:BOOL=OFF
        -DHDF5_BUILD_TOOLS:BOOL=OFF
        -DHDF5_BUILD_EXAMPLES:BOOL=OFF
        -DHDF5_BUILD_FORTRAN:BOOL=OFF
        -DHDF5_BUILD_JAVA:BOOL=OFF
        -DCMAKE_INSTALL_MESSAGE=NEVER #Avoid unnecessary output to console
        -DCMAKE_C_FLAGS=-w
        )

ExternalProject_Get_Property(HDF5 install_dir)
set(HDF5_LIB_DIR ${install_dir}/lib)
set(HDF5_INC_DIR ${install_dir}/include)
set(HDF5_SHARE_DIR ${install_dir}/share)
set(HDF5_CMAKE_DIR ${install_dir}/share/cmake)

file(WRITE  ${INSTALL_DIRECTORY}/hdf5/FindHDF5.cmake  "set(HDF5_LIB_DIR   ${HDF5_LIB_DIR})\n")
file(APPEND ${INSTALL_DIRECTORY}/hdf5/FindHDF5.cmake  "set(HDF5_INC_DIR   ${HDF5_INC_DIR})\n")
file(APPEND ${INSTALL_DIRECTORY}/hdf5/FindHDF5.cmake  "set(HDF5_SHARE_DIR ${HDF5_SHARE_DIR})\n")
file(APPEND ${INSTALL_DIRECTORY}/hdf5/FindHDF5.cmake  "set(HDF5_CMAKE_DIR ${HDF5_CMAKE_DIR})\n")
file(APPEND ${INSTALL_DIRECTORY}/hdf5/FindHDF5.cmake  "set(HDF5_INCLUDE_DIRS   ${HDF5_INC_DIR})\n") # For compatibility with find_package(HDF5...)






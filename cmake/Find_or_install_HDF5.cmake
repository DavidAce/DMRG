# The following script will find a native HDF5 installation on your system.
#
# If these packages are not found, they will be downloaded from a git repository and installed
# Either way, this script will define the following variables:
#
# HDF5_LIBRARIES              - Path to libhdf5*.a files
#
# HDF5_LIB_DIR                - Path do HDF5 library directory
# HDF5_INC_DIR                - Path to HDF5 include directory
# HDF5_BIN_DIR                - Path to HDF5 binary directory
#
# After execution, it is enough to target link libraries:
#   target_link_libraries(MyTarget ${HDF5_LIBRARIES})
# To use, simple include it in your CMakeLists.txt
#   include(Find_or_install_HDF5.cmake)
set(HDF5_USE_STATIC_LIBRARIES ON)
set(HDF5_FIND_DEBUG ON)
#set(HDF5_CXX_COMPILER_EXECUTABLE "h5c++")
find_package(HDF5 COMPONENTS CXX)

if (NOT HDF5_FOUND OR "${HDF5_LIBRARIES}" MATCHES "anaconda")
    find_package(PkgConfig)
    pkg_check_modules(HDF5 hdf5)
endif()


if (HDF5_FOUND)
    message("Found HDF5: ${HDF5_FOUND}")
    message("   In path: ${HDF5_LIBRARY_DIRS}")
    message("   LIBRARIES: ${HDF5_LIBRARIES}")
    message("   LDFLAGS:   ${HDF5_LDFLAGS}")
    message("   INCLUDE:   ${HDF5_INCLUDE_DIRS}")
endif()


if (HDF5_FOUND AND NOT "${HDF5_LIBRARIES}" MATCHES "anaconda")
        add_definitions    (${HDF5_DEFINITIONS})
        include_directories(${HDF5_INCLUDE_DIRS})
        message("FOUND NATIVE HDF5:   ${HDF5_LIBRARIES}")
else()
    message("DOWNLOADING HDF5...")
    execute_process(
            COMMAND ${CMAKE_COMMAND} -E make_directory tmp/hdf5
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts)
    execute_process(
            COMMAND ${CMAKE_COMMAND}
            -DINSTALL_DIRECTORY:PATH=${PROJECT_SOURCE_DIR}/libs
            -G ${CMAKE_GENERATOR} ../../hdf5
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/hdf5)
    execute_process(
            COMMAND ${CMAKE_COMMAND} --build .
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/hdf5
            RESULT_VARIABLE res_var)
    if(NOT "${res_var}" STREQUAL "0")
        message(FATAL_ERROR "HDF5 not found and failed to install: ${res_var}")
    endif()
    include(${PROJECT_SOURCE_DIR}/libs/hdf5/FindHDF5.cmake)
    set(HDF5_LIBRARIES "")
    get_libraries(${HDF5_LIB_DIR} hdf5  HDF5_LIBRARIES)
    include_directories(${HDF5_INC_DIR})
    message("SUCCESSFULLY INSTALLED HDF5:   ${HDF5_LIBRARIES}")
endif()



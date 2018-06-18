# The following script will find a native HDF5 installation on your system.
#
# If these packages are not found, they will be downloaded from a git repository and installed
# Either way, this script will define the following variables:
#
# HDF5_LIBRARIES              - Path to libhdf5*.a files
# HDF5_INCLUDE_DIR            - Path to libhdf5*.a files
#
# After execution, it is enough to target link libraries:
#   target_link_libraries(MyTarget ${HDF5_LIBRARIES})
# To use, simple include it in your CMakeLists.txt
#   include(Find_or_install_HDF5.cmake)

# Read more here: http://msarahan.github.io/cmake-and-hdf5-revisited/
message("SEARCHING FOR LIBRARY: HDF5")
if(EXISTS "${PROJECT_SOURCE_DIR}/libs/hdf5/FindHDF5.cmake")
    include(${PROJECT_SOURCE_DIR}/libs/hdf5/FindHDF5.cmake)
    set(HDF5_USE_STATIC_LIBRARIES ON)
    find_package(HDF5 COMPONENTS CXX HL static PATHS "${HDF5_CMAKE_DIR}" NO_DEFAULT_PATH NO_MODULE)
    if(HDF5_FOUND)
        list(APPEND HDF5_LIBRARIES ${HDF5_CXX_STATIC_LIBRARY} ${HDF5_HL_STATIC_LIBRARY} )
        list(REMOVE_DUPLICATES HDF5_LIBRARIES)
        target_include_directories(${PROJECT_NAME} PRIVATE ${HDF5_INCLUDE_DIR})
        target_link_libraries(${PROJECT_NAME} ${HDF5_LIBRARIES} -ldl)
        message(STATUS "FOUND PREVIOUSLY INSTALLED HDF5: ${HDF5_LIBRARIES} ")
    return()
    endif()
    message(WARNING "Found 'libs/hdf5/FindHDF5.cmake' in , but failed to configure. You may need a clean build.")
endif()



# Try finding in system
# Usually there is an executable present in system, "h5c++"
# Execute "find /usr -name h5c++ -print -quit" in terminal to find out where it is.
# -quit returns after first match
execute_process(
        COMMAND  find /usr -name libhdf5.* -exec dirname {} \; -quit
        OUTPUT_VARIABLE HDF5_ROOT
)

execute_process(
        COMMAND  find /usr -name h5c++ -print -quit
        OUTPUT_VARIABLE HDF5_CXX_COMPILER_EXECUTABLE
)

string(STRIP "${HDF5_CXX_COMPILER_EXECUTABLE}" HDF5_CXX_COMPILER_EXECUTABLE)
string(STRIP "${HDF5_ROOT}" HDF5_ROOT)

set(HDF5_USE_STATIC_LIBRARIES OFF)
set(HDF5_FIND_DEBUG OFF)
set(HDF5_DIR ${HDF5_ROOT})
message(STATUS "SEARCHING FOR HDF5 IN SYSTEM...")
find_package (HDF5 COMPONENTS CXX HL)
if (HDF5_FOUND AND NOT "${HDF5_LIBRARIES}" MATCHES "anaconda")
    add_definitions(${HDF5_DEFINITIONS})
    list(APPEND HDF5_LDFLAGS ${HDF5_CXX_LIBRARY_NAMES} ${HDF5_CXX_HL_LIBRARY_NAMES})
    list(APPEND HDF5_LIBRARIES ${HDF5_HL_LIBRARIES})
    list(REMOVE_DUPLICATES HDF5_LIBRARIES)

    message(STATUS "FOUND HDF5")
    message(STATUS "   In path: ${HDF5_INCLUDE_DIR}")
    message(STATUS "   HDF5 DEFINITIONS: ${HDF5_DEFINITIONS}")
    message(STATUS "   HDF5 LIBRARIES  : ${HDF5_LIBRARIES}")

    target_include_directories(${PROJECT_NAME} PRIVATE ${HDF5_INCLUDE_DIR})
    target_link_libraries(${PROJECT_NAME} ${HDF5_LIBRARIES} -ldl)
    return()
endif()






set(HDF5_FOUND 0)
set(HDF5_LIBRARIES "")
if (NOT HDF5_FOUND OR "${HDF5_LIBRARIES}" MATCHES "anaconda")
    message(STATUS "DOWNLOADING HDF5...")
    set(INSTALL_DIRECTORY ${PROJECT_SOURCE_DIR}/libs)

    execute_process(
            COMMAND ${CMAKE_COMMAND} -E make_directory tmp/hdf5
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts)
    execute_process(
            COMMAND ${CMAKE_COMMAND}
            -DINSTALL_DIRECTORY:PATH=${INSTALL_DIRECTORY}
            -G ${CMAKE_GENERATOR} ../../hdf5
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/hdf5)
    execute_process(
            OUTPUT_QUIET
            OUTPUT_FILE ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/hdf5/log_build.txt
            COMMAND ${CMAKE_COMMAND} --build .
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/hdf5
            RESULT_VARIABLE res_var)
    if(NOT "${res_var}" STREQUAL "0")
        message(FATAL_ERROR "HDF5 not found and failed to install: ${res_var}")
    endif()
    # Generate a file that can be detected if the library was installed successfully
    set(HDF5_CMAKE_DIR ${INSTALL_DIRECTORY}/hdf5/share/cmake)
    file(WRITE ${INSTALL_DIRECTORY}/hdf5/FindHDF5.cmake  "set(HDF5_CMAKE_DIR      ${HDF5_CMAKE_DIR})\n")


    #Include that file
    include(${INSTALL_DIRECTORY}/hdf5/FindHDF5.cmake)
    message("Searching in ${HDF5_CMAKE_DIR}")
    find_package(HDF5 COMPONENTS CXX HL PATHS "${HDF5_CMAKE_DIR}" NO_DEFAULT_PATH NO_MODULE REQUIRED static)
    if(HDF5_FOUND)
        list(APPEND HDF5_LIBRARIES ${HDF5_CXX_STATIC_LIBRARY} ${HDF5_HL_STATIC_LIBRARY} )
        message(STATUS "SUCCESSFULLY INSTALLED HDF5")
        message(STATUS "    In path: ${HDF5_INCLUDE_DIR}")
        message(STATUS "    HDF5 LIBRARIES  : ${HDF5_LIBRARIES}")
        message(STATUS "    BUILD LOG SAVED TO:   ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/hdf5/log_build.txt")
        target_include_directories(${PROJECT_NAME} PRIVATE ${HDF5_INCLUDE_DIR})
        target_link_libraries(${PROJECT_NAME} ${HDF5_LIBRARIES} -ldl)
        return()
    endif()
endif()


#
#    get_cmake_property(_variableNames VARIABLES)
#    foreach (_variableName ${_variableNames})
#        if("${_variableName}" MATCHES "HDF5")
#            message(STATUS "${_variableName}=${${_variableName}}")
#        endif()
#    endforeach()
# The following script will find a native BLAS and LAPACK installation on your system.
#
# If these packages are not found, they will be downloaded from a git repository and installed
# Either way, this script will define the following variables:
#
# LAPACK_LIBRARIES          - Path to liblapack.a
# BLAS_LIBRARIES            - Path to libblas.a, which comes included in the lapack installation.
#
# LAPACK_LIB                - Path do lapack library directory
# LAPACK_INC                - Path to lapack include directory
# LAPACK_BIN                - Path to lapack binary directory
#
# After execution, it is enough to target link libraries:
#   target_link_libraries(MyTarget ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
# To use, simple include it in your CMakeLists.txt
#   include(Find_or_install_LAPACK_and_BLAS.cmake)

#Check first if openblas is installed,
#set(BLA_VENDOR OpenBLAS)
#find_package(BLAS)

if(MKL_LIBRARIES)
    return()
endif()


message("SEARCHING FOR LIBRARY: LAPACK")
set(USE_OPTIMIZED_BLAS ON)
set(BLAS_FIND_QUIETLY ON)
if(EXISTS "${PROJECT_SOURCE_DIR}/libs/lapack/FindLapack.cmake" )
    include(cmake/FindGFortran.cmake)     ### For Fortran library
    include(${PROJECT_SOURCE_DIR}/libs/lapack/FindLapack.cmake)
    FILE(GLOB_RECURSE LAPACK_CMAKE_DIR "${LAPACK_CMAKE_DIR}/*lapack-config.cmake")
    get_filename_component(LAPACK_CMAKE_DIR "${LAPACK_CMAKE_DIR}" DIRECTORY)
    find_package(LAPACK PATHS "${LAPACK_CMAKE_DIR}" NO_DEFAULT_PATH NO_MODULE REQUIRED)
#    set(LAPACK_LIBRARIES ${LAPACK_lapack_LIBRARIES})
#    set(BLAS_LIBRARIES ${LAPACK_blas_LIBRARIES})
    get_property(LAPACK_LIBRARIES TARGET ${LAPACK_lapack_LIBRARIES} PROPERTY LOCATION)
    get_property(BLAS_LIBRARIES TARGET ${LAPACK_blas_LIBRARIES} PROPERTY LOCATION)
    if(LAPACK_lapack_LIBRARIES AND LAPACK_blas_LIBRARIES)
        message(STATUS "FOUND PREVIOUSLY INSTALLED LAPACK:   ${LAPACK_lapack_LIBRARIES} (Location -- ${LAPACK_LIBRARIES})")
        message(STATUS "FOUND PREVIOUSLY INSTALLED BLAS:     ${LAPACK_blas_LIBRARIES}   (Location -- ${BLAS_LIBRARIES})")
        target_link_libraries(${PROJECT_NAME} ${LAPACK_lapack_LIBRARIES})
        target_link_libraries(${PROJECT_NAME} ${LAPACK_blas_LIBRARIES})
        target_link_libraries(${PROJECT_NAME} ${GFORTRAN_LIB})
        return()
    endif()
    message(WARNING "Found 'libs/lapack/FindLAPACK.cmake' in , but failed to configure. You may need a clean build.")
endif()

message(STATUS "SEARCHING FOR LAPACK IN SYSTEM...")
find_package(LAPACK)
if(LAPACK_FOUND)
    message(STATUS "FOUND BLAS:     ${BLAS_LIBRARIES}")
    message(STATUS "FOUND LAPACK:   ${LAPACK_LIBRARIES}")
    target_link_libraries(${PROJECT_NAME} ${LAPACK_LIBRARIES})
    target_link_libraries(${PROJECT_NAME} ${BLAS_LIBRARIES})
else()
    message(STATUS "DOWNLOADING LAPACK...")
    set(INSTALL_DIRECTORY ../libs)
    include(cmake/FindGFortran.cmake)    ### For Fortran library
    unset(BLAS_LIBRARIES)
    unset(LAPACK_LIBRARIES)

    execute_process(
		COMMAND ${CMAKE_COMMAND} -E make_directory tmp/lapack
		WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts)
	execute_process(
		COMMAND ${CMAKE_COMMAND}
            -DINSTALL_DIRECTORY:PATH=${INSTALL_DIRECTORY}
            -G ${CMAKE_GENERATOR} ../../lapack
		WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/lapack )
	execute_process(
        OUTPUT_QUIET
        OUTPUT_FILE ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/lapack/log_build.txt
		COMMAND ${CMAKE_COMMAND} --build .
		WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/lapack
        RESULT_VARIABLE res_var)
    if(NOT "${res_var}" STREQUAL "0")
        message(FATAL_ERROR "LAPACK and BLAS not found and failed to install: ${res_var}")
    endif()
    # Generate a file that can be detected if the library was installed successfully
    set(LAPACK_CMAKE_DIR ${INSTALL_DIRECTORY}/lapack/lib/cmake)
    file(WRITE ${INSTALL_DIRECTORY}/lapack/FindLapack.cmake "set(LAPACK_CMAKE_DIR   ${LAPACK_CMAKE_DIR})\n")

    #Include that file
    include(${INSTALL_DIRECTORY}/lapack/FindLapack.cmake)
    FILE(GLOB_RECURSE LAPACK_CMAKE_DIR "${LAPACK_CMAKE_DIR}/*lapack-config.cmake")
    get_filename_component(LAPACK_CMAKE_DIR "${LAPACK_CMAKE_DIR}" DIRECTORY)
    find_package(LAPACK PATHS "${LAPACK_CMAKE_DIR}" NO_DEFAULT_PATH NO_MODULE REQUIRED)



#    set(LAPACK_LIBRARIES ${LAPACK_lapack_LIBRARIES})
#    set(BLAS_LIBRARIES ${LAPACK_blas_LIBRARIES})
    get_cmake_property(_variableNames VARIABLES)
    foreach (_variableName ${_variableNames})
        if("${_variableName}" MATCHES "BLAS" OR "${_variableName}" MATCHES "LAPACK" OR "${_variableName}" MATCHES "blas")
            message(STATUS "${_variableName}=${${_variableName}}")
        endif()
    endforeach()

    get_property(LAPACK_LIBRARIES TARGET ${LAPACK_lapack_LIBRARIES} PROPERTY LOCATION)
    get_property(BLAS_LIBRARIES TARGET ${LAPACK_blas_LIBRARIES} PROPERTY LOCATION)

    if(LAPACK_lapack_LIBRARIES AND LAPACK_blas_LIBRARIES)
        message(STATUS "SUCCESSFULLY INSTALLED LAPACK:   ${LAPACK_lapack_LIBRARIES}")
        message(STATUS "SUCCESSFULLY INSTALLED BLAS:     ${LAPACK_blas_LIBRARIES}")
        message(STATUS "BUILD LOG SAVED TO:   ${PROJECT_SOURCE_DIR}/cmake/download_scripts/tmp/lapack/log_build.txt")
        target_link_libraries(${PROJECT_NAME} ${LAPACK_lapack_LIBRARIES})
        target_link_libraries(${PROJECT_NAME} ${LAPACK_blas_LIBRARIES})
        target_link_libraries(${PROJECT_NAME} ${GFORTRAN_LIB})


#        get_property(blas_location TARGET blas PROPERTY LOCATION)
#        message (STATUS "blas_location == ${blas_location}")



    else()
        message(WARNING "LAPACK_LIBRARIES: ${LAPACK_LIBRARIES}")
        message(WARNING "BLAS_LIBRARIES: ${LAPACK_LIBRARIES}")
        message(FATAL_ERROR "LAPACK or BLAS could not be found.")
    endif()
endif()

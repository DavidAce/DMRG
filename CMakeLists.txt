cmake_minimum_required(VERSION 3.5)
project(Finite_DMRG_eigen)
set(CMAKE_VERBOSE_MAKEFILE OFF)

#######################
### Get git version ###
#######################
execute_process(
        COMMAND git rev-list HEAD --count
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_REVISION
        OUTPUT_STRIP_TRAILING_WHITESPACE
)
add_definitions("-DGIT_REVISION=${GIT_REVISION}")
message("Git revision: ${GIT_REVISION}")

# Require MPI for this project:
# find_package(MPI REQUIRED)
#include_directories(${MPI_INCLUDE_PATH})


# Require HDF5 for this project
########################################################
### IF HDF5 IS BUILT FROM SOURCE, USE THESE COMMANDS ###
########################################################
#set(HDF5_DIR "/usr/local/HDF5-1.10.1-Linux/HDF_Group/HDF5/1.10.1/share/cmake")
#find_package(HDF5 NAMES hdf5 COMPONENTS CXX CXX_HL NO_MODULE REQUIRED static)
#
#    link_libraries     (${HDF5_LIBRARIES})
#    link_libraries     (${HDF5_CXX_STATIC_LIBRARY})
#    link_libraries     (${HDF5_CXX_HL_STATIC_LIBRARY})
#    include_directories(${HDF5_INCLUDE_DIR})
#    add_definitions    (${HDF5_DEFINITIONS})
#
#
########################################################

set(HDF5_USE_STATIC_LIBRARIES ON)
find_package(HDF5 COMPONENTS CXX REQUIRED)
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()
if(HDF5_FOUND)
    message("Found HDF5:")
    message("   Version:               ${HDF5_VERSION}")
    message("   Included dirs:         ${HDF5_INCLUDE_DIRS}")
    include_directories(${HDF5_INCLUDE_DIRS})
    add_definitions    (${HDF5_DEFINITIONS})
    link_libraries     (${HDF5_LIBRARIES})
else()
    message("Please install HDF5")
    message(" > sudo apt install libhdf5-cpp-11 libhdf5-10 libhdf5-dev libhdf5-doc")
    message(FATAL_ERROR "HDF5 installation not found!")
endif()



#########################################
### Require GSL library               ###
#########################################
find_package(GSL REQUIRED)
link_libraries(${GSL_LIBRARIES})


# Require OpenMP
# find_package(OpenMP REQUIRED)

# include eigen3 and spectra
include_directories(eigen3)
include_directories(spectra/include)


#########################################
### Demand for c++17 compiler support ###
#########################################
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
message("Found compiler version: ${CMAKE_CXX_COMPILER_VERSION}. Standard: -std=cxx${CMAKE_CXX_STANDARD}")


######################################
### Change flags depending on host ###
######################################
site_name(this_host)
if ("${this_host}" MATCHES "triolith")
    message("Host: ${this_host} (Intel) ")
    set(CMAKE_CXX_COMPILER mpiicpc)
    set(CMAKE_CXX_FLAGS_DEBUG "  ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}     -g -O0 -ip -ipo -xavx -cxxlib=/software/apps/gcc/5.3.0/build01/ -gxx-name=g++ -Wall")
    set(CMAKE_CXX_FLAGS_RELEASE "  ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE} -g -O3 -ip -ipo -xavx -cxxlib=/software/apps/gcc/5.3.0/build01/ -gxx-name=g++ -DNDEBUG")
else()
    message("Host: ${this_host} (GNU)")
#    set(CMAKE_CXX_COMPILER mpicxx)
#    set(CMAKE_CXX_COMPILER g++)
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftime-report")

    if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-invalid-partial-specialization")
    endif()
    set(CMAKE_CXX_FLAGS_DEBUG " ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}     -Wall  -g -O0 -fopenmp -march=native -Wno-deprecated-declarations -Wno-ignored-attributes -Wno-int-in-bool-context")
    set(CMAKE_CXX_FLAGS_RELEASE " ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}        -g -O3 -fopenmp -march=native -Wno-deprecated-declarations -Wno-ignored-attributes -Wno-int-in-bool-context -DNDEBUG ")
endif()
message("Compiling with: ${CMAKE_CXX_COMPILER_ID}")


###########################################
###  Add all source files recursively   ###
###########################################
file(GLOB_RECURSE SOURCES "source/*.cpp")
file(GLOB_RECURSE HEADERS "source/*.h")

set (INCLUDE_DIRS "")
foreach (_headerFile ${HEADERS})
    get_filename_component(_dir ${_headerFile} PATH)
    list (APPEND INCLUDE_DIRS ${_dir})
endforeach()
list(REMOVE_DUPLICATES INCLUDE_DIRS)


include_directories(${INCLUDE_DIRS})
message("SOURCE FILES: ${SOURCES}")
add_executable(Finite_DMRG_eigen ${SOURCES} ${HEADERS} source/n_data_containers.h)

########################################################
###  target link libraries GO AFTER ADD_EXECUTABLE   ###
########################################################
target_link_libraries(Finite_DMRG_eigen -lstdc++fs) # For experimental::filesystem
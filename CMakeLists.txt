cmake_minimum_required(VERSION 3.5)
project(Finite_DMRG_eigen)
set(CMAKE_VERBOSE_MAKEFILE OFF)

# Require MPI for this project:
# find_package(MPI REQUIRED)
#include_directories(${MPI_INCLUDE_PATH})


# Require HDF5 for this project

########################################################
### IF HDF5 IS BUILT FROM SOURCE, USE THESE COMMANDS ###
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

#find_package(HDF5 COMPONENTS CXX HL REQUIRED)
find_package(HDF5 COMPONENTS CXX REQUIRED)
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()
if(HDF5_FOUND)
    message("Found HDF5:")
    message("   Version:               ${HDF5_VERSION}")
    message("   Included dirs:         ${HDF5_INCLUDE_DIRS}")
    message("   Required        libs:  ${HDF5_LIBRARIES}")
    message("   Required CXX    libs:  ${HDF5_CXX_LIBRARIES}")
    message("   Required CXX_HL libs:  ${HDF5_CXX_HL_LIBRARIES}")

    include_directories(${HDF5_INCLUDE_DIR})
    include_directories(${HDF5_INCLUDE_DIRS})
    add_definitions    (${HDF5_DEFINITIONS})
    link_libraries     (${HDF5_LIBRARIES})
    link_libraries     (${HDF5_HL_LIBRARIES})
    link_libraries     (${HDF5_CXX_LIBRARIES})
    link_libraries     (${HDF5_CXX_HL_LIBRARIES})
#    target_link_libraries(Finite_DMRG_eigen ${HDF5_CXX_LIBRARIES} ${HDF5_HL_LIBRARIES})
#    set(CMAKE_CXX_COMPILER ${HDF5_CXX_COMPILER_EXECUTABLE})
else()
    message("Please install HDF5")
    message(" > sudo apt install libhdf5-cpp-11 libhdf5-10 libhdf5-dev libhdf5-doc")
    message(FATAL_ERROR "HDF5 installation not found!")
endif()


# Require OpenMP
# find_package(OpenMP REQUIRED)

# include eigen3 and spectra
include_directories(eigen3)
include_directories(spectra/include)


##################################
# Check for c++17 compiler support
##################################
include(CheckCXXCompilerFlag)
set(ENABLE_CXXFLAGS_TO_CHECK
        -std=c++17
        -std=c++1z
        -std=c++14
        -std=c++1y
        -std=gnu++1y
        -std=gnu++14
        -std=c++11
        -std=gnu++11)

foreach(flag ${ENABLE_CXXFLAGS_TO_CHECK})
    string(REPLACE "-std=" "_" flag_var ${flag})
    string(REPLACE "+" "x" flag_var ${flag_var})
    check_cxx_compiler_flag("${flag}" COMPILER_HAS_CXX_FLAG${flag_var})
    if(COMPILER_HAS_CXX_FLAG${flag_var})
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}")
        message("Found compiler version: ${CMAKE_CXX_COMPILER_VERSION} ${flag}")
        break()
    endif()
endforeach()







###################################
#Change flags depending on host
###################################
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
add_executable(Finite_DMRG_eigen ${SOURCES} )
target_link_libraries(Finite_DMRG_eigen -lstdc++fs)
cmake_minimum_required(VERSION 3.5)
project(Finite_DMRG_eigen)
set(CMAKE_VERBOSE_MAKEFILE OFF)

# Require MPI for this project:
find_package(MPI REQUIRED)
#find_package(OpenMP REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
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
add_executable(Finite_DMRG_eigen ${SOURCES})


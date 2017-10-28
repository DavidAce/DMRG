cmake_minimum_required(VERSION 3.7)
include("hunter/HunterGate.cmake")
### Hunter snapshot that will be used ###
HunterGate(
        URL "https://github.com/ruslo/hunter/archive/v0.19.138.tar.gz"
        SHA1 "0b7c1f758843e73a8f0b3c71a58751377a617c22"
)
project(DMRG++)
set(CMAKE_VERBOSE_MAKEFILE OFF)

#######################
### Get git version ###
#######################

# Get the current working branch
execute_process(
        COMMAND git rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the latest abbreviated commit hash of the working branch
execute_process(
        COMMAND git log -1 --format=%h
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_COMMIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the revision count
execute_process(
        COMMAND git rev-list HEAD --count
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_REVISION
        OUTPUT_STRIP_TRAILING_WHITESPACE
)

configure_file(
${CMAKE_SOURCE_DIR}/gitversion/gitversion.h.in
${CMAKE_BINARY_DIR}/gitversion/gitversion.h
)

include_directories(${CMAKE_BINARY_DIR}/gitversion)


############################
### Import Eigen library ###
############################
hunter_add_package(Eigen)
find_package(Eigen3 CONFIG REQUIRED)


############################
### Import GSL library   ###
############################
hunter_add_package(GSL)
find_package(GSL CONFIG REQUIRED)

############################
### Import HDF5 library  ###
############################
hunter_add_package(hdf5)
find_package(ZLIB CONFIG REQUIRED)
find_package(szip CONFIG REQUIRED)
find_package(hdf5 CONFIG REQUIRED)


include_directories(spectra/include)


#########################################
### Demand for c++17 compiler support ###
#########################################
macro(use_cxx17)
  if (CMAKE_VERSION VERSION_LESS "3.8")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      set (CMAKE_CXX_FLAGS "-std=c++17 ${CMAKE_CXX_FLAGS}")
    endif ()
  else ()
    set (CMAKE_CXX_STANDARD 17)
  endif ()
endmacro(use_cxx17)
use_cxx17()
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
include_directories("source")
add_executable(DMRG++ ${SOURCES} ${HEADERS})



########################################################
###  target link libraries GO AFTER ADD_EXECUTABLE   ###
########################################################

##################################
### For #include <filesystem>  ###
##################################
target_link_libraries(DMRG++ -lstdc++fs) # For experimental::filesystem
target_link_libraries(DMRG++ Eigen3::Eigen)
target_link_libraries(DMRG++ GSL::gsl)
target_link_libraries(DMRG++ hdf5)


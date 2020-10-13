
message(STATUS "C compiler ${CMAKE_C_COMPILER}")
message(STATUS "FC compiler ${CMAKE_Fortran_COMPILER}")
message(STATUS "CXX compiler ${CMAKE_CXX_COMPILER}")

############################################################
### Set  the same microarchitecture for c++ and OpenBLAS ###
############################################################

if(NOT DMRG_MICROARCH)
    set(DMRG_MICROARCH "native")
endif()
if(DMRG_MICROARCH)
    if (${DMRG_MICROARCH} STREQUAL "zen")
        string(TOUPPER ${DMRG_MICROARCH} OPENBLAS_MARCH)
        set(CXX_MARCH znver1)
    elseif (${DMRG_MICROARCH} STREQUAL "native")
        set(OPENBLAS_MARCH HASWELL)
        set(CXX_MARCH native)
    else()
        string(TOUPPER ${DMRG_MICROARCH} OPENBLAS_MARCH)
        string(TOLOWER ${DMRG_MICROARCH} CXX_MARCH)
    endif()
endif()


######################################################################################
###                   Apply RELEASE/DEBUG compile flags                            ###
######################################################################################
# I have benchmarked the compiler flags below
#        -fstack-protector
#        -D_FORTIFY_SOURCE=2
#        -fno-omit-frame-pointer
#        -fno-strict-aliasing
# and found that there is NO significant performance difference in tensor contractions.
# These were the results from Tensorbench, best of 2 runs,
# variability between runs ~0.3 seconds.
# (all with -O3 -mfma -DNDEBUG -march=native -mtune=native):
# 19.3789 s: -fstack-protector
# 19.2768 s: -fstack-protector -fno-omit-frame-pointer -fno-strict-aliasing
# 19.6200 s: -fstack-protector -fno-omit-frame-pointer -fno-strict-aliasing -D_FORTIFY_SOURCE=2
# 19.3540 s: (none)
#
# In particular, -fno-strict-aliasing fixes a bug in tensor shuffle using -O3 -DNDEBUG
# in gcc 10.0.1. The bug causes the dimensions to blow up when assigning the shuffling
# operation to a new tensor, such as in:
#       Eigen::Tensor<Scalar, 3> mps_shuffled = mps.shuffle(Textra::array3{1, 0, 2});
# Since the dimensions become kind of random (and not near the limits of int/long),
# I believe gcc-10 is simply reordering something it shouldn't, thus assigning
# uninitialized values on these dimensions.
# The bug was hard to track down, and could be fixed by various other flags that
# hurt performance more, like -mno-avx, removing -DNDEBUG, or lowering -O3 to -O2.
######################################################################################

if(CXX_MARCH)
    message(STATUS "Using microarchitechture: ${CXX_MARCH}")
    list(APPEND CMAKE_CXX_FLAGS            -march=${CXX_MARCH} -mtune=${CXX_MARCH})
endif()

list(APPEND CMAKE_CXX_FLAGS                )
list(APPEND CMAKE_CXX_FLAGS_RELEASE  -g -O3 -fno-strict-aliasing -Wall -Wextra -Wpedantic)
list(APPEND CMAKE_CXX_FLAGS_DEBUG    -g -O0 -fno-strict-aliasing -Wall -Wextra -Wpedantic -fstack-protector -D_FORTIFY_SOURCE=2 -fno-omit-frame-pointer) #-D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC
list(APPEND CMAKE_CXX_FLAGS_RELWITHDEBINFO )
list(APPEND CMAKE_CXX_FLAGS_MINSIZEREL)

# For time tracing using ClangBuildAnalyzer
#if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
#    list(APPEND CMAKE_CXX_FLAGS -ftime-trace)
#endif()

string (REPLACE " " ";" CMAKE_CXX_FLAGS_LIST                "${CMAKE_CXX_FLAGS}")
string (REPLACE " " ";" CMAKE_CXX_FLAGS_RELEASE_LIST        "${CMAKE_CXX_FLAGS_RELEASE}")
string (REPLACE " " ";" CMAKE_CXX_FLAGS_DEBUG_LIST          "${CMAKE_CXX_FLAGS_DEBUG}")
string (REPLACE " " ";" CMAKE_CXX_FLAGS_RELWITHDEBINFO_LIST "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
string (REPLACE " " ";" CMAKE_CXX_FLAGS_MINSIZEREL_LIST     "${CMAKE_CXX_FLAGS_MINSIZEREL}")

list(REMOVE_DUPLICATES CMAKE_CXX_FLAGS_LIST)
list(REMOVE_DUPLICATES CMAKE_CXX_FLAGS_RELEASE_LIST)
list(REMOVE_DUPLICATES CMAKE_CXX_FLAGS_DEBUG_LIST)
list(REMOVE_DUPLICATES CMAKE_CXX_FLAGS_RELWITHDEBINFO_LIST)
list(REMOVE_DUPLICATES CMAKE_CXX_FLAGS_MINSIZEREL_LIST)

string (REPLACE ";" " " CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS_LIST}")
string (REPLACE ";" " " CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE_LIST}")
string (REPLACE ";" " " CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG_LIST}")
string (REPLACE ";" " " CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO_LIST}")
string (REPLACE ";" " " CMAKE_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_MINSIZEREL_LIST}")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND NOT CMAKE_EXE_LINKER_FLAGS MATCHES "fuse-ld=gold")
    set(CMAKE_EXE_LINKER_FLAGS "-fuse-ld=gold ${CMAKE_EXE_LINKER_FLAGS}")
endif()


###########################################
###  Transmit compiler to dependencies  ###
###########################################

if(NOT DEFINED ENV{CC} AND DEFINED CMAKE_C_COMPILER)
    set(ENV{CC} ${CMAKE_C_COMPILER})
endif()
if(NOT DEFINED ENV{CXX} AND DEFINED CMAKE_CXX_COMPILER)
    set(ENV{CXX} ${CMAKE_CXX_COMPILER})
endif()
if(NOT DEFINED ENV{FC} AND DEFINED CMAKE_Fortran_COMPILER)
    set(ENV{FC} ${CMAKE_Fortran_COMPILER})
endif()


###############################
# Settings for shared builds

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already (but later on when installing)
# Note: Since DMRG++ is often run from the build folder we want to keep the build-folder RPATH in the executable.
#       Therefore itt makes sense to keep this setting "FALSE" here but "TRUE" for dependencies that are
#       installed with in "fetch" mode with externalproject_add
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)



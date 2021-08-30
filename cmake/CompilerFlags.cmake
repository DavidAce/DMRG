
message(DEBUG "C compiler ${CMAKE_C_COMPILER}")
message(DEBUG "FC compiler ${CMAKE_Fortran_COMPILER}")
message(DEBUG "CXX compiler ${CMAKE_CXX_COMPILER}")

#####################################################
### Set the  microarchitecture for OpenBLAS       ###
#####################################################
cmake_host_system_information(RESULT _host_name QUERY HOSTNAME)
set(OPENBLAS_TARGET HASWELL)
set(OPENBLAS_DYNAMIC_ARCH ON)
if($ENV{CI} OR $ENV{GITHUB_ACTIONS} OR DMRG_MICROARCH MATCHES "generic|Generic|GENERIC")
    set(MARCH -march=x86-64)
    set(MTUNE -mtune=generic)
    set(OPENBLAS_TARGET GENERIC)
    set(OPENBLAS_DYNAMIC_ARCH OFF)
elseif(DEFINED DMRG_MICROARCH)
    set(MARCH -march=${DMRG_MICROARCH})
    set(MTUNE -mtune=${DMRG_MICROARCH})
else()
    set(MARCH -march=haswell)
    set(MTUNE -mtune=native)
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

#######################################
### Write compile commands to file  ###
#######################################
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


message(DEBUG "Using ${MARCH} ${MTUNE}")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${MARCH} ${MTUNE}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -fno-strict-aliasing -Wall -Wextra -Wpedantic")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-strict-aliasing -Wall -Wextra -Wpedantic -fstack-protector -D_FORTIFY_SOURCE=2 -fno-omit-frame-pointer") #-D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC
#set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND NOT CMAKE_EXE_LINKER_FLAGS MATCHES "fuse-ld=gold")
    set(CMAKE_EXE_LINKER_FLAGS "-fuse-ld=gold -Wl,--disable-new-dtags")
endif()

# Set these variables so that the same flags are used for building dependencies
set(CMAKE_CXX_FLAGS_INIT                 "${CMAKE_CXX_FLAGS_INIT} ${CMAKE_CXX_FLAGS}" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE_INIT         "${CMAKE_CXX_FLAGS_RELEASE_INIT} ${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_DEBUG_INIT           "${CMAKE_CXX_FLAGS_DEBUG_INIT} ${CMAKE_CXX_FLAGS_DEBUG}" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT  "${CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT} ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}" CACHE STRING "" FORCE)
set(CMAKE_EXE_LINKER_FLAGS_INIT          "${CMAKE_EXE_LINKER_FLAGS_INIT} ${CMAKE_EXE_LINKER_FLAGS}" CACHE STRING "" FORCE)





# For time tracing using ClangBuildAnalyzer
#if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
#    list(APPEND CMAKE_CXX_FLAGS -ftime-trace)
#endif()

###############################
# Settings for shared builds
###############################

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already (but later on when installing)
# Note: Since DMRG++ is often run from the build folder we want to keep the build-folder RPATH in the executable.
#       Therefore it makes sense to keep this setting "FALSE" here but "TRUE" for dependencies that are
#       installed with in "cmake" mode with externalproject_add
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)



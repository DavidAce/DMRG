cmake_minimum_required(VERSION 3.10)
project(DMRG++ CXX )
#enable_language()

set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake-modules)


################################################################
### Set options for compilation                              ###
################################################################
option(BUILD_SHARED_LIBS  "Sets shared linking on/off"                          OFF)
option(BUILD_EXAMPLES     "Build examples"                                      OFF)
option(ENABLE_TESTS       "Enable testing"                                      ON )
option(DOWNLOAD_MISSING   "Download and install missing dependencies"           ON )
option(ENABLE_CCACHE      "Use ccache to speed up compilation"                  OFF)
option(ENABLE_COTIRE      "Use cotire to speed up compilation"                  OFF)
option(ENABLE_OPENMP      "Use OpenMP"                                          OFF)
option(USE_MKL            "Use Intel MKL libraries"                             ON )
option(ENABLE_LTO         "Enable link time optimization"                       OFF)

#################################################################
### Set default policies if CMake is new enough               ###
#################################################################
if (CMAKE_VERSION VERSION_LESS 3.12)
    message(STATUS "Not setting policies")
else()
    cmake_policy(SET CMP0074 NEW)
    cmake_policy(SET CMP0075 NEW)
    # Use cmake_policy CMP0067 to let try_compile honor current CMAKE_CXX_STANDARD setting.
    # See here https://stackoverflow.com/questions/47213356/cmake-using-corrext-c-standard-when-checking-for-header-files
    cmake_policy(SET CMP0067 NEW)
endif()


if(USE_MKL AND NOT MKLROOT)
    if(DEFINED $ENV{MKLROOT})
        set(MKLROOT $ENV{MKLROOT})
    else()
    set(MKLROOT  /opt/intel/mkl)
    message(STATUS "MKLROOT is not defined. Setting default: ${MKLROOT}")
    endif()
endif()


################################################################
### Get git version number                                   ###
### Generates a header gitversion/gitversion.h               ###
### Include it using #include <gitversion.h>                 ###
### Gives a namespace GIT:: with several git version numbers.###
################################################################
include(cmake-modules/gitversion.cmake)



##################################################################
### Get operating system properties into a variable            ###
##################################################################
cmake_host_system_information(RESULT _host_name QUERY HOSTNAME)

set(GET_OS_INFO_CMD lsb_release -a)
if(${CMAKE_HOST_APPLE})
    set(GET_OS_INFO_CMD "sw_vers")
endif()
execute_process(COMMAND ${GET_OS_INFO_CMD}
        OUTPUT_VARIABLE OS_PROPERTIES
        OUTPUT_STRIP_TRAILING_WHITESPACE)
message("========= DETECTED OS =========")
message("Hostname: " ${_host_name})
message("${OS_PROPERTIES}")
message("===============================")




###################
### Setup flags ###
###################
message("Compiling with: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")

if(NOT MARCH)
    set(MARCH "native")
endif()

if(MARCH)
    if (${MARCH} STREQUAL "zen")
        string(TOUPPER ZEN OPENBLAS_MARCH)
        string(TOLOWER znver1 CXX_MARCH)
    elseif (${MARCH} STREQUAL "native")
        set(OPENBLAS_MARCH "haswell")
        string(TOLOWER ${MARCH} CXX_MARCH)
    else()
        string(TOUPPER ${MARCH} OPENBLAS_MARCH)
        string(TOLOWER ${MARCH} CXX_MARCH)
    endif()
endif()


if(ENABLE_CCACHE)
    find_program(CCACHE_FOUND ccache)
    if(CCACHE_FOUND)
        set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
        set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
        set(CMAKE_CXX_COMPILER_LAUNCHER ccache)
    endif(CCACHE_FOUND)
endif()

#################################
### Microarch dependent flags ###
#################################
# Check supported options here http://man7.org/linux/man-pages/man1/gcc.1.html
if(${MARCH} STREQUAL "haswell" OR ${MARCH} STREQUAL "skylake" OR ${MARCH} STREQUAL "zen" )
    #list(APPEND RELEASE_OPTIONS -mmmx -msse -msse2 -msse3 -mssse3 -msse4.2 -msse4.1 -mpopcnt -maes -mpclmul -mavx -mavx2 -mmovbe -mrdrnd -mf16c -mfsgsbase -mfma  -mbmi -mbmi2 -fstack-protector-strong -Wformat -Wformat-security)
else()
    # Set default to whatever sandybridge can manage
    #list(APPEND RELEASE_OPTIONS  -mmmx -msse -msse2 -msse3 -mssse3 -msse4.2 -msse4.1  -mpopcnt  -maes -mpclmul -mavx -fstack-protector-strong -Wformat -Wformat-security)
endif()











###########################################
###  Add all source files recursively   ###
###########################################
unset(SOURCES CACHE)
unset(HEADERS CACHE)
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "source/*.cpp")
file(GLOB_RECURSE HEADERS CONFIGURE_DEPENDS "source/*.h")
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})
target_include_directories(${PROJECT_NAME} PUBLIC source)

################################################
###  Force cmake to find .a library suffixes ###
################################################
if(BUILD_SHARED_LIBS)
    set(CUSTOM_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CUSTOM_SUFFIX} ${CMAKE_FIND_LIBRARY_SUFFIXES})
else()
    set(CUSTOM_SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CUSTOM_SUFFIX} ${CMAKE_FIND_LIBRARY_SUFFIXES})
    target_link_libraries  (${PROJECT_NAME} PRIVATE -static)                                             ### Static linkage
endif()





###############################
###  Require c++17 support  ###
###############################
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)


#######################################################
###  Check  #include<filesystem> or                 ###
###  or     #include<experimental/filesystem>       ###
#######################################################
include(${PROJECT_SOURCE_DIR}/cmake-modules/CheckCXXFilesystem.cmake)
CheckCXXFilesystem()

#######################################################
###  Check  #include<experimental/type_traits>      ###
#######################################################
include(${PROJECT_SOURCE_DIR}/cmake-modules/CheckTypeTraits.cmake)
CheckTypeTraits()

#######################################################
###  Check  #include<optional> or                   ###
###  or     #include<experimental/optional>         ###
#######################################################
include(${PROJECT_SOURCE_DIR}/cmake-modules/CheckCXXOptional.cmake)
CheckCXXOptional()


if (NOT TARGET Threads::Threads)
    ### Adapt pthread for static/dynamic linking
    ### This is to fix a bug with threads described here:
    ### https://stackoverflow.com/questions/35116327/when-g-static-link-pthread-cause-segmentation-fault-why
    set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
    set(THREADS_PREFER_PTHREAD_FLAG FALSE)
    find_package(Threads)
    include(cmake-modules/PrintTargetProperties.cmake)
    print_target_properties(Threads::Threads)
    if(TARGET Threads::Threads AND NOT BUILD_SHARED_LIBS)
        set_target_properties(Threads::Threads PROPERTIES INTERFACE_LINK_LIBRARIES "-Wl,--whole-archive -lpthread -Wl,--no-whole-archive")
    endif()
endif()




###########################################
###  Apply RELEASE/DEBUG compile flags  ###
###########################################
add_library(project-settings INTERFACE) # A helper target
message(STATUS "Using microarchitechture: ${CXX_MARCH}")
target_compile_options(project-settings INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-std=c++17 -g -Wall -Wpedantic -Wextra  -march=${CXX_MARCH} -mtune=${CXX_MARCH}>)
target_compile_options(project-settings INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:RELEASE>>:          -O0 >)
target_compile_options(project-settings INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:DEBUG>>:            -O0 -g3 -fstack-protector -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -D_FORTIFY_SOURCE=2>)
target_compile_options(project-settings INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:RELWITHDEBINFO>>:   -O1 -g3 -fstack-protector -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -D_FORTIFY_SOURCE=2>)
target_compile_options(project-settings INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:PROFILE>>:          -O3 -DNDEBUG  -lprofiler -g  -ftime-report>)


#######################################
### Compiler-dependent linker flags ###
#######################################
target_link_libraries(project-settings INTERFACE  -lstdc++fs)
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" )
    target_link_libraries(project-settings INTERFACE  -stdlib=libstdc++ )
endif()

target_compile_features(${PROJECT_NAME}  INTERFACE cxx_std_17)
set_target_properties(${PROJECT_NAME} PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO
        )

###  Optional OpenMP support
###  Note that Clang has some  trouble with static openmp and that
###  and that static openmp is not recommended. This tries to enable
###  static openmp anyway because I find it useful. Installing
###  libiomp5 might help for shared linking.
if(ENABLE_OPENMP)
    include(cmake-modules/FindPackageOpenMP.cmake)
    find_package_openmp()
    if(TARGET OpenMP)
        target_link_libraries(project-settings INTERFACE OpenMP)
    else()
        target_compile_options(project-settings INTERFACE -Wno-unknown-pragmas)
    endif()
endif()

##################################################################
### Check if we can enable Link Time Optimization              ###
##################################################################
if(ENABLE_LTO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT supported OUTPUT error)
    if(supported)
        message(STATUS "IPO / LTO supported, enabling")
        set_target_properties(${PROJECT_NAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION ON)
    else()
        message(STATUS "IPO / LTO not supported: <${error}>, will not enable")
    endif()
endif()


##################################################################
### Setup dependency paths                                     ###
##################################################################
include(GNUInstallDirs)
string(TOLOWER ${CMAKE_BUILD_TYPE} cmake_build_type)
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    message(STATUS "Setting default install prefix -- CMAKE_INSTALL_PREFIX --> ${PROJECT_SOURCE_DIR}/libs-${cmake_build_type}")
    set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/libs-${cmake_build_type}" CACHE PATH "..." FORCE)
endif()
set(EXTERNAL_INSTALL_DIR       ${CMAKE_INSTALL_PREFIX})
set(EXTERNAL_BUILD_DIR         ${CMAKE_BINARY_DIR}/external-deps)
set(EXTERNAL_CONFIG_DIR        ${CMAKE_BINARY_DIR}/external-deps)





##########################################
###  Find or fetch all dependencies   ###
##########################################



include(cmake-modules/FindGFortran.cmake)

if(USE_MKL)
    include(cmake-modules/Find_dont_install_INTELMKL.cmake)    # MKL - Intel's math Kernel Library, use the BLAS implementation in Eigen and Arpack. Includes lapack.
endif()

include(cmake-modules/Fetch_OpenBLAS.cmake)                 # If MKL is not found openblas will be used instead. Includes lapack.
include(cmake-modules/FindLAPACKE.cmake)                    # If MKL is not found openblas will be used instead. Includes lapack.
include(cmake-modules/Fetch_Arpack-ng.cmake)                # Iterative Eigenvalue solver for a few eigenvalues/eigenvectors using Arnoldi method.
include(cmake-modules/Fetch_Arpackpp.cmake)                 # LC++ frontend for arpack-ng
include(cmake-modules/Fetch_Eigen3.cmake)                   # Eigen3 numerical library
include(cmake-modules/Fetch_h5pp.cmake)                     # GSL for computing some numerical integrals.
include(cmake-modules/Fetch_Gflags.cmake)
include(cmake-modules/Fetch_Glog.cmake)
include(cmake-modules/Fetch_Ceres.cmake)

target_link_libraries(${PROJECT_NAME} PRIVATE arpack++)
target_link_libraries(${PROJECT_NAME} PRIVATE lapacke)
target_link_libraries(${PROJECT_NAME} PRIVATE Eigen3::Eigen )
target_link_libraries(${PROJECT_NAME} PRIVATE h5pp::h5pp h5pp::deps h5pp::flags)
target_link_libraries(${PROJECT_NAME} PRIVATE ceres)
target_link_libraries(${PROJECT_NAME} PRIVATE project-settings) # This needs to go last in case OpenMP is enabled
target_link_libraries(${PROJECT_NAME} PRIVATE Threads::Threads)

##################################################################
### Print summary of CMake configuration                       ###
##################################################################
message("=============================================================================")
message("===               DMRG++: SUMMARY OF CMAKE CONFIGURATION                  ===")
message("=============================================================================")
include(cmake-modules/PrintTargetInfo.cmake)
print_target_info(OpenMP)
print_target_info(blas)
print_target_info(lapack)
print_target_info(lapacke)
print_target_info(arpack)
print_target_info(arpack++)
print_target_info(Eigen3::Eigen)
print_target_info(h5pp::h5pp)
print_target_info(h5pp::deps)
print_target_info(h5pp::flags)
print_target_info(spdlog::spdlog)
print_target_info(hdf5::hdf5)
print_target_info(gflags)
print_target_info(glog::glog)
print_target_info(ceres)
print_target_info(project-settings)
print_target_info(Threads::Threads)
message("=============================================================================")
message("")


###########################################
### Use precompiled headers with cotire ###
###########################################
if (ENABLE_COTIRE)
    get_target_property(LAPACKE_DIRS   lapacke         INTERFACE_INCLUDE_DIRECTORIES)
    get_target_property(ARPACK++_DIRS  arpack++        INTERFACE_INCLUDE_DIRECTORIES)
    set_target_properties(${PROJECT_NAME} PROPERTIES
            COTIRE_PREFIX_HEADER_INCLUDE_PATH "${PROJECT_SOURCE_DIR}/source"
            COTIRE_PREFIX_HEADER_IGNORE_PATH
            "${LAPACKE_DIRS};${ARPACK++_DIRS}"
            )
    include(cmake-modules/cotire.cmake)
    cotire(${PROJECT_NAME})

endif()



###################
### INFORMATION ###
###################
message("")
message("=============================================================================")
message("")
message("      To clean cmake files, run './build.sh -c'")
message("      To clean cmake files downloaded libraries, run './build.sh -l '")
message("      To launch this program,  './run.sh'")
message("")
message("=============================================================================")
message("")



#############################
## UNIT TEST
#############################
if(ENABLE_TESTS)
    enable_testing()
    add_subdirectory(tests/arpack++_simple)
    add_subdirectory(tests/arpack++_mps   )
    add_subdirectory(tests/file_io        )
    add_subdirectory(tests/eigsolver      )
    add_subdirectory(tests/svdwrapper     )
    add_subdirectory(tests/eigen_openmp   )
    add_subdirectory(tests/ceres_pthread  )

    add_custom_target(all-tests ALL
            DEPENDS test-arpack++_simple
            DEPENDS test-arpack++_mps
            DEPENDS test-file_io
            DEPENDS test-eigsolver
            DEPENDS test-eigen_openmp
            DEPENDS test-svdwrapper
            DEPENDS test-ceres_pthread
            )

    add_custom_command(TARGET all-tests
            POST_BUILD
            COMMENT "Running Tests"
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            COMMAND ${CMAKE_CTEST_COMMAND} -C $<CONFIGURATION> --output-on-failures)
endif()






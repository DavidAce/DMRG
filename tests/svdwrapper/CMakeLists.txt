
cmake_minimum_required(VERSION 3.10)
#project(svdwrapper)
add_executable(TEST-svdwrapper svdwrapper.cpp
        ${PROJECT_SOURCE_DIR}/source/math/class_svd_wrapper.h
        ${PROJECT_SOURCE_DIR}/source/math/class_svd_wrapper.cpp
        ${PROJECT_SOURCE_DIR}/source/math/class_svd_wrapper_lapacke.cpp)

set(CMAKE_VERBOSE_MAKEFILE ON)

target_include_directories(TEST-svdwrapper PUBLIC ${PROJECT_SOURCE_DIR}/source)
target_include_directories(TEST-svdwrapper PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(TEST-svdwrapper PRIVATE TEST_DIR="${PROJECT_SOURCE_DIR}/tests/svdwrapper")

################################################
###  Force cmake to find .a library suffixes ###
################################################
if(NOT BUILD_SHARED_LIBS)
    target_link_libraries  (TEST-svdwrapper PRIVATE -static)                                             ### Static linkage
endif()
target_link_libraries(TEST-svdwrapper PRIVATE Eigen3::Eigen lapacke h5pp::h5pp  hdf5   spdlog::spdlog)
target_link_libraries(TEST-svdwrapper PRIVATE ${COMPILER_LFLAGS})



set_target_properties  (TEST-svdwrapper PROPERTIES CXX_STANDARD_REQUIRED 17)
target_compile_features(TEST-svdwrapper PRIVATE cxx_std_17)
target_compile_options (TEST-svdwrapper PRIVATE ${COMMON_OPTIONS})                                   ### Common options
target_compile_options (TEST-svdwrapper PRIVATE "$<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>")               ### Debug build options
target_compile_options (TEST-svdwrapper PRIVATE "$<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>")           ### Release build options
target_compile_options (TEST-svdwrapper PRIVATE "$<$<CONFIG:RELWITHDEBINFO>:${RELWITHDEBINFO_OPTIONS}>")    ### Release build with debug info options
add_test(NAME svdwrapper_test COMMAND TEST-svdwrapper)

add_dependencies(TEST-svdwrapper h5pp::h5pp hdf5 Eigen3::Eigen spdlog::spdlog gfortran )
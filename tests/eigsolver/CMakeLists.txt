
cmake_minimum_required(VERSION 3.10)
#project(eigsolver)
add_executable(TEST-eigsolver eigsolver.cpp
        ${PROJECT_SOURCE_DIR}/source/general/class_eigsolver.cpp
        ${PROJECT_SOURCE_DIR}/source/general/class_tic_toc.cpp)


target_include_directories(TEST-eigsolver PUBLIC ${PROJECT_SOURCE_DIR}/source)
target_include_directories(TEST-eigsolver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(TEST-eigsolver PRIVATE TEST_MATRIX_DIR="${PROJECT_SOURCE_DIR}/tests/eigsolver")

################################################
###  Force cmake to find .a library suffixes ###
################################################
if(NOT BUILD_SHARED_LIBS)
    target_link_libraries  (TEST-eigsolver PRIVATE -static)                                             ### Static linkage
endif()

target_link_libraries(TEST-eigsolver PRIVATE -lstdc++fs)
target_link_libraries(TEST-eigsolver PRIVATE -flto)
target_link_libraries(TEST-eigsolver PRIVATE Eigen3::Eigen h5pp::h5pp  hdf5   spdlog::spdlog)
target_link_libraries(TEST-eigsolver PRIVATE gfortran)

set_target_properties  (TEST-eigsolver PROPERTIES CXX_STANDARD_REQUIRED 17)
target_compile_features(TEST-eigsolver PRIVATE cxx_std_17)
target_compile_options (TEST-eigsolver PRIVATE ${COMMON_OPTIONS})                                   ### Common options
target_compile_options (TEST-eigsolver PRIVATE "$<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>")               ### Debug build options
target_compile_options (TEST-eigsolver PRIVATE "$<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>")           ### Release build options
add_test(NAME eigsolver_test COMMAND TEST-eigsolver)

add_dependencies(TEST-eigsolver h5pp::h5pp hdf5 Eigen3::Eigen spdlog::spdlog )
cmake_minimum_required(VERSION 3.10)
#project(file_io)

add_executable(TEST-file_io file_io.cpp
        ../../source/general/nmspc_tensor_extra.h
        )
set(CMAKE_VERBOSE_MAKEFILE ON)

target_include_directories(TEST-file_io PUBLIC ${PROJECT_SOURCE_DIR}/source)
target_include_directories(TEST-file_io PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
################################################
###  Force cmake to find .a library suffixes ###
################################################
if(NOT BUILD_SHARED_LIBS)
    target_link_libraries  (TEST-file_io PRIVATE -static)                                             ### Static linkage
endif()


target_link_libraries(TEST-file_io PRIVATE Eigen3::Eigen)
target_link_libraries(TEST-file_io PRIVATE hdf5::hdf5   )
target_link_libraries(TEST-file_io PRIVATE spdlog::spdlog)
target_link_libraries(TEST-file_io PRIVATE h5pp::h5pp)
target_link_libraries(TEST-file_io PRIVATE ${COMPILER_LFLAGS})


set_target_properties  (TEST-file_io PROPERTIES CXX_STANDARD_REQUIRED 17)
target_compile_features(TEST-file_io PRIVATE cxx_std_17)
target_compile_options (TEST-file_io PRIVATE "${COMMON_OPTIONS}")                                        ### Common options
target_compile_options (TEST-file_io PRIVATE "$<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>")                      ### Debug build options
target_compile_options (TEST-file_io PRIVATE "$<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>")                  ### Release build options
target_compile_options (TEST-file_io PRIVATE "$<$<CONFIG:RELWITHDEBINFO>:${RELWITHDEBINFO_OPTIONS}>")    ### Release build with debug info options

add_test(NAME file_io_test COMMAND TEST-file_io)
#add_dependencies(TEST-file_io Eigen3 hdf5::hdf5 hdf5::hdf5_hl hdf5::hdf5_cpp hdf5::hdf5_hl_cpp)
add_dependencies(TEST-file_io Eigen3::Eigen hdf5::hdf5 spdlog::spdlog h5pp::h5pp)
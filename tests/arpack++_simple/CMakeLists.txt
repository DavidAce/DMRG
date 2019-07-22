
cmake_minimum_required(VERSION 3.10)
#project(arpack++_simple)

add_executable(TEST-arpack++_simple arpack++_simple.cpp)
set(CMAKE_VERBOSE_MAKEFILE ON)

################################################
###  Force cmake to find .a library suffixes ###
################################################
if(NOT BUILD_SHARED_LIBS)
    target_link_libraries  (TEST-arpack++_simple PRIVATE -static)                                             ### Static linkage
endif()

target_link_libraries(TEST-arpack++_simple PRIVATE  arpack++)
target_link_libraries(TEST-arpack++_simple PRIVATE ${COMPILER_LFLAGS})


add_dependencies(TEST-arpack++_simple arpack++)


set_target_properties  (TEST-arpack++_simple PROPERTIES CXX_STANDARD_REQUIRED 17)
target_compile_features(TEST-arpack++_simple PRIVATE cxx_std_17)
target_compile_options (TEST-arpack++_simple PRIVATE ${COMMON_OPTIONS})                                          ### Common options
target_compile_options (TEST-arpack++_simple PRIVATE "$<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>")                      ### Debug build options
target_compile_options (TEST-arpack++_simple PRIVATE "$<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>")                  ### Release build options
target_compile_options (TEST-arpack++_simple PRIVATE "$<$<CONFIG:RELWITHDEBINFO>:${RELWITHDEBINFO_OPTIONS}>")    ### Release build with debug info options

add_test(NAME arpack++_simple_test COMMAND TEST-arpack++_simple)



cmake_minimum_required(VERSION 3.10)
#project(arpack++_mps)

enable_language(Fortran)

add_executable(TEST-arpack++_mps arpack++_mps.cpp
                            ../../source/general/class_eigsolver.cpp
                            ../../source/general/arpack_extra/arpackpp_solver.cpp
                            ../../source/general/arpack_extra/matrix_product_stl.cpp
                            ../../source/general/class_tic_toc.cpp)
target_include_directories(TEST-arpack++_mps PUBLIC ${PROJECT_SOURCE_DIR}/source)
target_include_directories(TEST-arpack++_mps PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

################################################
###  Force cmake to find .a library suffixes ###
################################################
if(NOT BUILD_SHARED_LIBS)
    target_link_libraries  (TEST-arpack++_mps PRIVATE -static)                                             ### Static linkage
endif()


target_link_libraries(TEST-arpack++_mps PRIVATE  arpack++ Eigen3::Eigen spdlog::spdlog)
#######################################
### Compiler-dependent linker flags ###
#######################################
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU" )
    target_link_libraries  (${PROJECT_NAME} PRIVATE -lstdc++fs)                                          ### For experimental::filesystem
elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    target_link_libraries  (${PROJECT_NAME} PRIVATE -stdlib=libc++)                                      ### For experimental::filesystem
endif()


set_target_properties  (TEST-arpack++_mps PROPERTIES CXX_STANDARD_REQUIRED 17)
target_compile_features(TEST-arpack++_mps PRIVATE cxx_std_17)
target_compile_options (TEST-arpack++_mps PRIVATE "${COMMON_OPTIONS}")                                 ### Common options
target_compile_options (TEST-arpack++_mps PRIVATE "$<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>")               ### Debug build options
target_compile_options (TEST-arpack++_mps PRIVATE "$<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>")           ### Release build options

add_test(NAME arpack++_mps_test COMMAND TEST-arpack++_mps)
add_dependencies(TEST-arpack++_mps arpack++ Eigen3::Eigen spdlog::spdlog)
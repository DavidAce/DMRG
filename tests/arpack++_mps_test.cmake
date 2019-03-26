enable_language(Fortran)

add_executable(arpack++_mps_test_target tests/arpack++_mps_test.cpp
                                        source/general/class_eigsolver.cpp
                                        source/general/arpack_extra/arpackpp_solver.cpp
                                        source/general/arpack_extra/matrix_product_stl.cpp
                                        source/general/class_tic_toc.cpp)
set_target_properties(arpack++_mps_test_target PROPERTIES OUTPUT_NAME  arpack++_mps_test)
target_include_directories(arpack++_mps_test_target PUBLIC source)


################################################
###  Force cmake to find .a library suffixes ###
################################################
if(NOT BUILD_SHARED_LIBS)
    target_link_libraries  (arpack++_mps_test_target PRIVATE -static)                                             ### Static linkage
endif()


target_link_libraries(arpack++_mps_test_target PRIVATE  arpack++ Eigen3::Eigen spdlog::spdlog)
target_link_libraries(arpack++_mps_test_target PRIVATE -lstdc++fs)
target_link_libraries(arpack++_mps_test_target PRIVATE -flto)

set_target_properties  (arpack++_mps_test_target PROPERTIES CXX_STANDARD_REQUIRED 17)
target_compile_features(arpack++_mps_test_target PRIVATE cxx_std_17)
target_compile_options (arpack++_mps_test_target PRIVATE "${COMMON_OPTIONS}")                                 ### Common options
target_compile_options (arpack++_mps_test_target PRIVATE "$<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>")               ### Debug build options
target_compile_options (arpack++_mps_test_target PRIVATE "$<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>")           ### Release build options

add_test(NAME arpack++_mps_test COMMAND arpack++_mps_test_target)
add_dependencies(arpack++_mps_test_target arpack++ Eigen3::Eigen spdlog::spdlog)
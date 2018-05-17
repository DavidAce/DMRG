add_executable(arpack++_mps_test_target tests/arpack++_mps_test.cpp source/general/class_arpack_eigsolver.cpp source/general/class_arpack_custom_products.cpp)
set_target_properties(arpack++_mps_test_target PROPERTIES OUTPUT_NAME  arpack++_mps_test_object)
target_link_libraries(arpack++_mps_test_target EIGEN3 arpack arpackpp  ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${EXTRA_LDLAGS})
target_include_directories(arpack++_mps_test_target PRIVATE ${arpack++_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} source)
set_target_properties  (arpack++_mps_test_target PROPERTIES CXX_STANDARD_REQUIRED 17)
target_compile_features(arpack++_mps_test_target PRIVATE cxx_std_17)
target_compile_options (arpack++_mps_test_target PRIVATE ${COMMON_OPTIONS} ${EIGEN3_COMPILER_FLAGS})                                   ### Common options
target_compile_options (arpack++_mps_test_target PRIVATE "$<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>")               ### Debug build options
target_compile_options (arpack++_mps_test_target PRIVATE "$<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>")           ### Release build options
add_test(NAME arpack++_mps_test COMMAND arpack++_mps_test_target)


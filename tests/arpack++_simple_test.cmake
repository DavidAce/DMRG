add_executable(arpack++_simple_test_target tests/arpack++_simple_test.cpp)
set_target_properties(arpack++_simple_test_target PROPERTIES OUTPUT_NAME  arpack++_simple_test_object)
target_link_libraries(arpack++_simple_test_target arpack arpackpp ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${EXTRA_LDLAGS})
target_include_directories(arpack++_simple_test_target PRIVATE ${arpack++_INCLUDE_DIR})
add_test(NAME arpack++_simple_test COMMAND arpack++_simple_test_target)


enable_language(Fortran)

get_target_property(BLAS_LIBRARIES        blas INTERFACE_LINK_LIBRARIES    )
get_target_property(BLAS_INCLUDE_DIRS     blas INTERFACE_INCLUDE_DIRECTORY )
get_target_property(BLAS_LINK_FLAGS       blas INTERFACE_LINK_FLAGS        )
get_target_property(BLAS_COMPILE_FLAGS    blas INTERFACE_COMPILE_OPTIONS   )
#
get_target_property(LAPACK_LIBRARIES        lapack INTERFACE_LINK_LIBRARIES    )
get_target_property(LAPACK_INCLUDE_DIRS     lapack INTERFACE_INCLUDE_DIRECTORY )
get_target_property(LAPACK_LINK_FLAGS       lapack INTERFACE_LINK_FLAGS        )
get_target_property(LAPACK_COMPILE_FLAGS    lapack INTERFACE_COMPILE_OPTIONS   )

add_executable(arpack++_simple_test_target tests/arpack++_simple_test.cpp)

set_target_properties(arpack++_simple_test_target PROPERTIES OUTPUT_NAME  arpack++_simple_test_object)
target_link_libraries(arpack++_simple_test_target PRIVATE arpack arpackpp blas lapack ${BLAS_LINK_FLAGS} ${LAPACK_LINK_FLAGS})
target_include_directories(arpack++_simple_test_target PRIVATE ${arpack++_INCLUDE_DIR} ${BLAS_INCLUDE_DIRS} ${LAPACK_INCLUDE_DIRS})
target_compile_options(arpack++_simple_test_target PRIVATE ${BLAS_COMPILE_FLAGS} ${LAPACK_COMPILE_FLAGS})
set_target_properties  (arpack++_simple_test_target PROPERTIES CXX_STANDARD_REQUIRED 17)
target_compile_features(arpack++_simple_test_target PRIVATE cxx_std_17)
target_compile_options (arpack++_simple_test_target PRIVATE ${COMMON_OPTIONS})                                   ### Common options
target_compile_options (arpack++_simple_test_target PRIVATE "$<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>")               ### Debug build options
target_compile_options (arpack++_simple_test_target PRIVATE "$<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>")           ### Release build options

add_test(NAME arpack++_simple_test COMMAND arpack++_simple_test_target)


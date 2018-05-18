enable_language(Fortran)

add_executable(arpack++_mps_test_target tests/arpack++_mps_test.cpp source/general/class_arpack_eigsolver.cpp source/general/class_arpack_custom_products.cpp)
set_target_properties(arpack++_mps_test_target PROPERTIES OUTPUT_NAME  arpack++_mps_test_object)


target_link_libraries(arpack++_mps_test_target PRIVATE -v arpack arpack++ blas lapack EIGEN3)
target_compile_options(arpack++_mps_test_target
        PRIVATE -v
        $<TARGET_PROPERTY:arpack,INTERFACE_COMPILE_OPTIONS>
        $<TARGET_PROPERTY:arpack++,INTERFACE_COMPILE_OPTIONS>
        $<TARGET_PROPERTY:blas,INTERFACE_COMPILE_OPTIONS>
        $<TARGET_PROPERTY:lapack,INTERFACE_COMPILE_OPTIONS>
        $<TARGET_PROPERTY:EIGEN3,INTERFACE_COMPILE_OPTIONS>
        )

target_include_directories(arpack++_mps_test_target
        PRIVATE
        source
        $<TARGET_PROPERTY:arpack,INTERFACE_INCLUDE_DIRECTORY>
        $<TARGET_PROPERTY:arpack++,INTERFACE_INCLUDE_DIRECTORY>
        $<TARGET_PROPERTY:blas,INTERFACE_INCLUDE_DIRECTORY>
        $<TARGET_PROPERTY:lapack,INTERFACE_INCLUDE_DIRECTORY>
        $<TARGET_PROPERTY:EIGEN3,INTERFACE_INCLUDE_DIRECTORY>
        )

set_target_properties  (arpack++_mps_test_target PROPERTIES CXX_STANDARD_REQUIRED 17)
target_compile_features(arpack++_mps_test_target PRIVATE cxx_std_17)
target_compile_options (arpack++_mps_test_target PRIVATE ${COMMON_OPTIONS})                                   ### Common options
target_compile_options (arpack++_mps_test_target PRIVATE "$<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>")               ### Debug build options
target_compile_options (arpack++_mps_test_target PRIVATE "$<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>")           ### Release build options

add_test(NAME arpack++_mps_test COMMAND arpack++_mps_test_target)
add_dependencies(arpack++_mps_test_target blas lapack arpack arpack++ EIGEN3)
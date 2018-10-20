add_executable(arpack++_simple_test_target tests/arpack++_simple_test.cpp)
set_target_properties(arpack++_simple_test_target PROPERTIES OUTPUT_NAME  arpack++_simple_test)
target_link_libraries(arpack++_simple_test_target  PRIVATE  arpack++ arpack blas lapack gfortran)
#target_link_libraries(arpack++_simple_test_target PRIVATE arpack++)
#target_link_libraries(arpack++_simple_test_target PRIVATE arpack)
#target_link_libraries(arpack++_simple_test_target PRIVATE blas)
#target_link_libraries(arpack++_simple_test_target PRIVATE lapack)
#target_link_libraries(arpack++_simple_test_target PRIVATE ${QUADMATH_LIB})
#target_link_libraries(arpack++_simple_test_target PRIVATE ${GFORTRAN_LIB})
target_link_libraries(arpack++_simple_test_target PRIVATE -lpthread)
#target_link_libraries(arpack++_simple_test_target PRIVATE -lomp)
target_link_libraries(arpack++_simple_test_target PRIVATE -lstdc++fs)
target_link_libraries(arpack++_simple_test_target PRIVATE -flto)
target_compile_options(arpack++_simple_test_target
        PRIVATE
        $<TARGET_PROPERTY:arpack,INTERFACE_COMPILE_OPTIONS>
        $<TARGET_PROPERTY:arpack++,INTERFACE_COMPILE_OPTIONS>
        $<TARGET_PROPERTY:blas,INTERFACE_COMPILE_OPTIONS>
        $<TARGET_PROPERTY:lapack,INTERFACE_COMPILE_OPTIONS>
        )

target_include_directories(arpack++_simple_test_target
        PRIVATE
        $<TARGET_PROPERTY:arpack,INTERFACE_INCLUDE_DIRECTORY>
        $<TARGET_PROPERTY:arpack++,INTERFACE_INCLUDE_DIRECTORY>
        $<TARGET_PROPERTY:blas,INTERFACE_INCLUDE_DIRECTORY>
        $<TARGET_PROPERTY:lapack,INTERFACE_INCLUDE_DIRECTORY>
        )

set_target_properties  (arpack++_simple_test_target PROPERTIES CXX_STANDARD_REQUIRED 17)
target_compile_features(arpack++_simple_test_target PRIVATE cxx_std_17)
target_compile_options (arpack++_simple_test_target PRIVATE ${COMMON_OPTIONS})                                   ### Common options
target_compile_options (arpack++_simple_test_target PRIVATE -m64)                                                ### 64 bit binary
target_compile_options (arpack++_simple_test_target PRIVATE "$<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>")               ### Debug build options
target_compile_options (arpack++_simple_test_target PRIVATE "$<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>")           ### Release build options
add_test(NAME arpack++_simple_test COMMAND arpack++_simple_test_target)

add_dependencies(arpack++_simple_test_target  arpack arpack++ blas lapack gfortran)
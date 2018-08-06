enable_language(Fortran)

add_executable(arpack++_mps_test_target tests/arpack++_mps_test.cpp
                                        source/general/class_eigsolver_arpack.cpp
                                        source/general/class_eigsolver_arpack_custom_products.cpp
                                        source/general/class_tic_toc.cpp)
set_target_properties(arpack++_mps_test_target PROPERTIES OUTPUT_NAME  arpack++_mps_test)
target_link_libraries(arpack++_mps_test_target PRIVATE arpack arpack++ blas lapack EIGEN3 -lstdc++fs -flto)

set_target_properties  (arpack++_mps_test_target PROPERTIES CXX_STANDARD_REQUIRED 17)
target_compile_features(arpack++_mps_test_target PRIVATE cxx_std_17)
target_compile_options (arpack++_mps_test_target PRIVATE "${COMMON_OPTIONS};-Wno-unused-but-set-variable")                                   ### Common options
target_compile_options (arpack++_mps_test_target PRIVATE "$<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>")               ### Debug build options
target_compile_options (arpack++_mps_test_target PRIVATE "$<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>")           ### Release build options


target_compile_options(arpack++_mps_test_target
        PRIVATE
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

add_test(NAME arpack++_mps_test COMMAND arpack++_mps_test_target)
add_dependencies(${PROJECT_NAME} arpack++_mps_test_target blas lapack arpack arpack++ EIGEN3 elemental)
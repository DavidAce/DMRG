add_executable(sanity_test_target tests/sanity_test.cpp)


set_target_properties(sanity_test_target PROPERTIES OUTPUT_NAME  sanity_test)
target_include_directories(sanity_test_target PUBLIC source)
if(NOT BUILD_SHARED_LIBS)
    target_link_libraries  (sanity_test_target PRIVATE -static)                                             ### Static linkage
endif()
target_link_libraries(sanity_test_target PRIVATE -flto)
target_link_libraries(sanity_test_target PRIVATE blas spdlog::spdlog Eigen3::Eigen)

set_target_properties  (sanity_test_target PROPERTIES CXX_STANDARD_REQUIRED 17)
target_compile_features(sanity_test_target PRIVATE cxx_std_17)
target_compile_options (sanity_test_target PRIVATE ${COMMON_OPTIONS})                                   ### Common options
target_compile_options (sanity_test_target PRIVATE "$<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>")               ### Debug build options
target_compile_options (sanity_test_target PRIVATE "$<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>")           ### Release build options
add_test(NAME sanity_test COMMAND sanity_test_target)

add_dependencies(sanity_test_target  blas )
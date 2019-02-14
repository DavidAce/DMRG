
add_executable(hdf5_test_target tests/hdf5_test.cpp
        source/IO/class_hdf5_file.cpp
        source/IO/class_hdf5_file.h
        source/general/nmspc_tensor_extra.h
        )
set_target_properties(hdf5_test_target PROPERTIES OUTPUT_NAME  hdf5_test)
target_include_directories(hdf5_test_target PUBLIC source)

target_link_libraries(hdf5_test_target PRIVATE eigen3)
target_link_libraries(hdf5_test_target PRIVATE hdf5      )
target_link_libraries(hdf5_test_target PRIVATE spdlog    )
#target_link_libraries(hdf5_test_target PRIVATE hdf5::hdf5_hl    )
#target_link_libraries(hdf5_test_target PRIVATE hdf5::hdf5_cpp   )
#target_link_libraries(hdf5_test_target PRIVATE hdf5::hdf5_hl_cpp)
target_link_libraries(hdf5_test_target PRIVATE -lstdc++fs)
target_link_libraries(hdf5_test_target PRIVATE -flto)

set_target_properties  (hdf5_test_target PROPERTIES CXX_STANDARD_REQUIRED 17)
target_compile_features(hdf5_test_target PRIVATE cxx_std_17)
target_compile_options (hdf5_test_target PRIVATE "${COMMON_OPTIONS}")                                 ### Common options
target_compile_options (hdf5_test_target PRIVATE "$<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>")               ### Debug build options
target_compile_options (hdf5_test_target PRIVATE "$<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>")           ### Release build options

add_test(NAME hdf5_test COMMAND hdf5_test_target)
#add_dependencies(hdf5_test_target eigen3 hdf5::hdf5 hdf5::hdf5_hl hdf5::hdf5_cpp hdf5::hdf5_hl_cpp)
add_dependencies(hdf5_test_target eigen3 hdf5 )
add_executable(arpack++_simple_test_target tests/arpack++_simple_test.cpp)
set_target_properties(arpack++_simple_test_target PROPERTIES OUTPUT_NAME  arpack++_simple_test)

################################################
###  Force cmake to find .a library suffixes ###
################################################
if(STATIC_BUILD)
    set(CUSTOM_SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CUSTOM_SUFFIX} ${CMAKE_FIND_LIBRARY_SUFFIXES})
    target_link_libraries  (arpack++_simple_test_target PRIVATE -static)                  ### Static linkage
else()
    set(CUSTOM_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CUSTOM_SUFFIX} ${CMAKE_FIND_LIBRARY_SUFFIXES})
endif()

target_link_libraries  (arpack++_simple_test_target PRIVATE -Wl,--no-allow-shlib-undefined )
target_link_libraries  (arpack++_simple_test_target PRIVATE -Wl,--no-as-needed )
target_link_libraries  (arpack++_simple_test_target PRIVATE -Wl,--no-undefined )

target_link_libraries(arpack++_simple_test_target PRIVATE -lstdc++fs)
target_link_libraries(arpack++_simple_test_target PRIVATE -flto)
target_link_libraries(arpack++_simple_test_target PRIVATE  arpack++)


set_target_properties  (arpack++_simple_test_target PROPERTIES CXX_STANDARD_REQUIRED 17)
target_compile_features(arpack++_simple_test_target PRIVATE cxx_std_17)
target_compile_options (arpack++_simple_test_target PRIVATE ${COMMON_OPTIONS})                                   ### Common options
target_compile_options (arpack++_simple_test_target PRIVATE "$<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>")               ### Debug build options
target_compile_options (arpack++_simple_test_target PRIVATE "$<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>")           ### Release build options
add_test(NAME arpack++_simple_test COMMAND arpack++_simple_test_target)

add_dependencies(arpack++_simple_test_target arpack++)
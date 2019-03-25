add_executable(lapacke_test_target tests/lapacke_test.cpp
        source/general/class_eigsolver.h
        source/general/class_eigsolver.cpp
        source/general/arpack_extra/arpackpp_solver.cpp
        source/general/arpack_extra/matrix_product_stl.cpp
        source/general/class_tic_toc.cpp)


set_target_properties(lapacke_test_target PROPERTIES OUTPUT_NAME  lapacke_test)
target_include_directories(lapacke_test_target PUBLIC source)

################################################
###  Force cmake to find .a library suffixes ###
################################################
if(STATIC_BUILD)
    set(CUSTOM_SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CUSTOM_SUFFIX} ${CMAKE_FIND_LIBRARY_SUFFIXES})
    target_link_libraries  (lapacke_test_target PRIVATE -static)                  ### Static linkage
else()
    set(CUSTOM_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CUSTOM_SUFFIX} ${CMAKE_FIND_LIBRARY_SUFFIXES})
endif()

target_link_libraries  (lapacke_test_target PRIVATE -Wl,--no-allow-shlib-undefined )
target_link_libraries  (lapacke_test_target PRIVATE -Wl,--no-as-needed )
target_link_libraries  (lapacke_test_target PRIVATE -Wl,--no-undefined )

target_link_libraries(lapacke_test_target PRIVATE -lstdc++fs)
target_link_libraries(lapacke_test_target PRIVATE -flto)
target_link_libraries(lapacke_test_target PRIVATE  arpack++  h5pp::h5pp hdf5 Eigen3::Eigen spdlog::spdlog LBFGSpp)


set_target_properties  (lapacke_test_target PROPERTIES CXX_STANDARD_REQUIRED 17)
target_compile_features(lapacke_test_target PRIVATE cxx_std_17)
target_compile_options (lapacke_test_target PRIVATE ${COMMON_OPTIONS})                                   ### Common options
target_compile_options (lapacke_test_target PRIVATE "$<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>")               ### Debug build options
target_compile_options (lapacke_test_target PRIVATE "$<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>")           ### Release build options
add_test(NAME lapacke_test COMMAND lapacke_test_target)

add_dependencies(lapacke_test_target arpack++  h5pp::h5pp hdf5 Eigen3::Eigen spdlog::spdlog LBFGSpp)
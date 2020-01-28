
cmake_minimum_required(VERSION 3.10)
project(test-arpack++_mps)

add_executable(test-arpack++_mps arpack++_mps.cpp
                            ../../source/math/class_eigsolver.cpp
                            ../../source/math/arpack_extra/arpackpp_solver.cpp
                            ../../source/math/arpack_extra/matrix_product_hamiltonian.cpp
                            ../../source/math/arpack_extra/matrix_product_dense.cpp
                            ../../source/math/arpack_extra/matrix_product_sparse.cpp
                            ../../source/math/arpack_extra/matrix_product_stl.cpp
                            ../../source/general/class_tic_toc.cpp

        )
target_include_directories(test-arpack++_mps PRIVATE ${CMAKE_SOURCE_DIR}/source)
target_include_directories(test-arpack++_mps PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(test-arpack++_mps PRIVATE project-settings )
add_test(NAME test-arpack++_mps COMMAND test-arpack++_mps)


cmake_minimum_required(VERSION 3.10)
project(test-arpack++_mps)

add_executable(test-arpack++_mps arpack++_mps.cpp
                            ../../source/math/class_eigsolver.cpp
                            ../../source/math/arpack_extra/arpackpp_solver.cpp
                            ../../source/math/arpack_extra/matrix_product_hamiltonian.cpp
                            ../../source/math/arpack_extra/matrix_product_hamiltonian_sq.cpp
                            ../../source/math/arpack_extra/matrix_product_dense.cpp
                            ../../source/math/arpack_extra/matrix_product_sparse.cpp
                            ../../source/math/arpack_extra/matrix_product_stl.cpp
                            ../../source/general/class_tic_toc.cpp

        )
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/source)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(${PROJECT_NAME} PRIVATE project-settings )
set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS "/ignore:4099")
add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})

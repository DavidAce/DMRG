add_library(dmrg-math-eig-primme OBJECT)
target_sources(dmrg-math-eig-primme PRIVATE solver_primme.cpp)
target_link_libraries(dmrg-math-eig-primme PUBLIC deps flags)
target_link_precompiled_headers(dmrg-math-eig-primme)

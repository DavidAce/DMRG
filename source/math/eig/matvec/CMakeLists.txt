add_library(dmrg-math-eig-matvec OBJECT)
target_sources(dmrg-math-eig-matvec PUBLIC matvec_dense.cpp matvec_sparse.cpp matvec_mpo.cpp)
target_link_libraries(dmrg-math-eig-matvec PUBLIC deps flags)
target_link_precompiled_headers(dmrg-math-eig-matvec)


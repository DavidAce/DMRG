add_library(dmrg-math-svd OBJECT)
target_sources(dmrg-math-svd PRIVATE svd.cpp svd_eigen.cpp svd_lapacke.cpp svd_rsvd.cpp)
target_link_libraries(dmrg-math-svd PUBLIC deps flags)
target_link_precompiled_headers(dmrg-math-svd)

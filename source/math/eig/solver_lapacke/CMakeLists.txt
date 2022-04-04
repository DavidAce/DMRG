# This one is special, since it absolutely can't have BLAS/LAPACK headers during compile for some reason
add_library(dmrg-math-eig-lapacke OBJECT)
target_sources(dmrg-math-eig-lapacke PRIVATE
        solver_dgeev.cpp
        solver_dsyevd.cpp
        solver_dsyevx.cpp
        solver_dsyevr.cpp
        solver_zgeev.cpp
        solver_zheevd.cpp)
target_link_libraries(dmrg-math-eig-lapacke PUBLIC deps flags)
target_link_precompiled_headers(dmrg-math-eig-lapacke)

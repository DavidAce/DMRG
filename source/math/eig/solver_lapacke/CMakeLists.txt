# This one is special, since it absolutely can't have BLAS/LAPACK headers during compile for some reason
add_library(dmrg-math-eig-lapacke OBJECT)
target_sources(dmrg-math-eig-lapacke PRIVATE
        solver_dgeev.cpp
        solver_dsyevd.cpp
        solver_zgeev.cpp
        solver_zheevd.cpp)
target_link_libraries(dmrg-math-eig-lapacke PUBLIC dmrg-deps dmrg-flags)
if (DMRG_ENABLE_PCH)
    target_precompile_headers(dmrg-math-eig-lapacke REUSE_FROM dmrg-pch)
endif ()

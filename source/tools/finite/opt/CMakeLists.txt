add_library(dmrg-tools-finite-opt OBJECT)
target_sources(dmrg-tools-finite-opt PRIVATE
        ceres_base.cpp
        ceres_direct.cpp
        ceres_direct_functor.cpp
        ceres_subspace.cpp
        ceres_subspace_find_full.cpp
        ceres_subspace_find_part.cpp
        ceres_subspace_filter.cpp
        ceres_subspace_functor.cpp
        ceres_subspace_hamiltonian.cpp
        ground_state.cpp
        krylov_optimization.cpp
        report.cpp
        opt.cpp
        opt_mps.cpp
        opt_meta.cpp
        )

target_link_libraries(dmrg-tools-finite-opt PUBLIC dmrg-deps dmrg-flags)
if (DMRG_ENABLE_PCH)
    target_precompile_headers(dmrg-tools-finite-opt REUSE_FROM dmrg-pch)
endif ()
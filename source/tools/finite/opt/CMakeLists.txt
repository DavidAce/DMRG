add_library(dmrg-tools-finite-opt OBJECT)
target_sources(dmrg-tools-finite-opt PRIVATE
        ceres_base.cpp
        ceres_direct.cpp
        ceres_direct_functor.cpp
        ceres_subspace.cpp
        ceres_subspace_diag_full.cpp
        ceres_subspace_diag_part.cpp
        ceres_subspace_filter.cpp
        ceres_subspace_functor.cpp
        ceres_subspace_hamiltonian.cpp
        ground_state.cpp
        krylov_energy.cpp
        krylov_variance.cpp
        krylov_solution.cpp
        report.cpp
        opt.cpp
        opt_mps.cpp
        )

target_link_libraries(dmrg-tools-finite-opt PUBLIC dmrg-deps dmrg-flags)
if (DMRG_ENABLE_PCH)
    target_precompile_headers(dmrg-tools-finite-opt REUSE_FROM dmrg-pch)
endif ()
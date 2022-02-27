add_library(dmrg-tools-finite-opt OBJECT)
target_sources(dmrg-tools-finite-opt PRIVATE
        bfgs_base_functor.cpp
        bfgs_optimize_variance.cpp
        bfgs_variance_functor.cpp
        bfgs_subspace_functor.cpp
        bfgs_simps_functor.cpp
        bfgs_callback.cpp
        eigs_optimize_subspace.cpp
        eigs_extract_results.cpp
        eigs_filter_subspace.cpp
        eigs_optimize_energy.cpp
        eigs_optimize_variance.cpp
        eigs_optimize_overlap.cpp
        eigs_optimize_subspace.cpp
        eigs_find_subspace.cpp
        report.cpp
        opt.cpp
        opt_mps.cpp
        opt_meta.cpp
        )

target_link_libraries(dmrg-tools-finite-opt PUBLIC deps flags)
target_link_precompiled_headers(dmrg-tools-finite-opt)

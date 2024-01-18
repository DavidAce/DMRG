add_source_target(dmrg-tools-finite-opt
        TARGET_SOURCES
            find_subspace.cpp
            filter_subspace.cpp
            optimize_overlap.cpp
            optimize_energy.cpp
            optimize_energy_eig.cpp
            optimize_variance.cpp
            optimize_variance_subspace.cpp
            optimize_variance_eig.cpp
            opt.cpp
            opt_mps.cpp
            opt_meta.cpp
            opt_extract_results.cpp
            report.cpp
        OBJECT_LINK_LIBRARIES
            dmrg-deps dmrg-flags
        INTERFACE_LINK_LIBRARIES
            dmrg-tid
            dmrg-config
            dmrg-math-eig
            dmrg-algorithms-status
            dmrg-tensors-finite
            dmrg-tools-finite-measure
            dmrg-tools-common-log
            dmrg-tools-common-contraction
        )


add_source_target(dmrg-tools-infinite-measure
        TARGET_SOURCES
        measure.cpp
        measure_ham.cpp
        measure_mom.cpp
        OBJECT_LINK_LIBRARIES
        dmrg-deps dmrg-flags
        INTERFACE_LINK_LIBRARIES
        dmrg-tid
        dmrg-math-tenx
        dmrg-algorithms-status
        dmrg-tensors-infinite
        dmrg-tools-common-contraction
        dmrg-tools-common-log
        )


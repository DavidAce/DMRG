add_source_target(dmrg-tools-finite-ops
        TARGET_SOURCES
            ops.cpp
        OBJECT_LINK_LIBRARIES
            dmrg-deps dmrg-flags
        INTERFACE_LINK_LIBRARIES
            dmrg-tid
            dmrg-qm-spin
            dmrg-qm-mpo
            dmrg-tensors-state-finite
            dmrg-tensors-mpo
            dmrg-math-tenx
            dmrg-tools-finite-mps
            dmrg-tools-finite-measure
            dmrg-tools-common-log
            dmrg-tools-common-contraction
        )

add_source_target(dmrg-tools-finite-mpo
        TARGET_SOURCES
            mpo.cpp mpo_inv.cpp
        OBJECT_LINK_LIBRARIES
            dmrg-deps dmrg-flags
        INTERFACE_LINK_LIBRARIES
            dmrg-tid
            dmrg-math-svd
            dmrg-math-tenx
            dmrg-tensors-mpo
            dmrg-tools-common-log
            dmrg-tools-common-contraction
        )
add_source_target(dmrg-tools-finite-print
        TARGET_SOURCES
            print.cpp
        OBJECT_LINK_LIBRARIES
            dmrg-deps dmrg-flags
        INTERFACE_LINK_LIBRARIES
            dmrg-tensors-finite
            dmrg-tools-common-log
        )


add_subdirectory(h5)
add_subdirectory(ed)
add_subdirectory(measure)
add_subdirectory(multisite)
add_subdirectory(env)
add_subdirectory(mps)
add_subdirectory(opt)

add_library(dmrg-tools-finite INTERFACE)
target_link_libraries(dmrg-tools-finite INTERFACE
        dmrg-tools-finite-h5
        dmrg-tools-finite-ed
        dmrg-tools-finite-measure
        dmrg-tools-finite-mps
        dmrg-tools-finite-env
        dmrg-tools-finite-multisite
        dmrg-tools-finite-opt
        dmrg-tools-finite-mpo
        dmrg-tools-finite-ops
        dmrg-tools-finite-print
        )
target_link_libraries(dmrg-tools-finite INTERFACE dmrg-deps dmrg-flags)

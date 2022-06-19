add_source_target(dmrg-algorithms-status
        TARGET_SOURCES
        AlgorithmStatus.cpp
        OBJECT_LINK_LIBRARIES
        deps flags
        )


add_source_target(dmrg-algorithms-base
        TARGET_SOURCES
        AlgorithmBase.cpp
        OBJECT_LINK_LIBRARIES
        deps flags
        INTERFACE_LINK_LIBRARIES
        dmrg-config
        dmrg-algorithms-status
        dmrg-tools-common-h5
        dmrg-tools-common-log
        )

add_source_target(dmrg-algorithms-finite
        TARGET_SOURCES
        AlgorithmFinite.cpp
        OBJECT_LINK_LIBRARIES
        deps flags
        INTERFACE_LINK_LIBRARIES
        dmrg-tid
        dmrg-config
        dmrg-debug
        dmrg-qm-spin
        dmrg-tensors-state-finite
        dmrg-tensors-model-finite
        dmrg-tensors-edges-finite
        dmrg-algorithms-status
        dmrg-algorithms-base
        dmrg-tools-finite-env
        dmrg-tools-finite-h5
        dmrg-tools-finite-measure
        dmrg-tools-finite-ops
        dmrg-tools-finite-print
        dmrg-tools-common-log
        )

add_source_target(dmrg-algorithms-infinite
        TARGET_SOURCES
        AlgorithmInfinite.cpp
        OBJECT_LINK_LIBRARIES
        deps flags
        INTERFACE_LINK_LIBRARIES
        dmrg-tid
        dmrg-config
        dmrg-debug
        dmrg-algorithms-status
        dmrg-algorithms-base
        dmrg-tensors-state-finite
        dmrg-tensors-model-finite
        dmrg-tools-infinite-h5
        dmrg-tools-infinite-measure
        dmrg-tools-common-log
        dmrg-tools-common-h5

        )
add_source_target(dmrg-algorithms-fdmrg
        TARGET_SOURCES
        fdmrg.cpp
        OBJECT_LINK_LIBRARIES
        deps flags
        INTERFACE_LINK_LIBRARIES
        dmrg-algorithms-finite
        dmrg-tools-finite-opt
        dmrg-tools-common-h5
        dmrg-tools-common-prof

        )
add_source_target(dmrg-algorithms-flbit
        TARGET_SOURCES
        flbit.cpp
        OBJECT_LINK_LIBRARIES
        deps flags
        INTERFACE_LINK_LIBRARIES
        dmrg-algorithms-finite
        dmrg-math-tenx
        dmrg-qm-lbit
        dmrg-tools-finite-opt
        dmrg-tools-finite-mps
        dmrg-tools-common-h5
        dmrg-tools-common-prof
        )

add_source_target(dmrg-algorithms-xdmrg
        TARGET_SOURCES
        xdmrg.cpp
        OBJECT_LINK_LIBRARIES
        deps flags
        INTERFACE_LINK_LIBRARIES
        dmrg-algorithms-finite
        dmrg-qm-time
        dmrg-math-rnd
        dmrg-algorithms-fdmrg
        dmrg-tools-finite-opt
        )

add_source_target(dmrg-algorithms-idmrg
        TARGET_SOURCES
        idmrg.cpp
        OBJECT_LINK_LIBRARIES
        deps flags
        INTERFACE_LINK_LIBRARIES
        dmrg-algorithms-infinite
        dmrg-tensors-edges-infinite
        dmrg-tools-infinite-opt
        )

add_source_target(dmrg-algorithms-itebd
        TARGET_SOURCES
        itebd.cpp
        OBJECT_LINK_LIBRARIES
        deps flags
        INTERFACE_LINK_LIBRARIES
        dmrg-qm-time
        dmrg-algorithms-infinite
        dmrg-tools-infinite-opt

        )

add_source_target(dmrg-algorithms-launcher
        TARGET_SOURCES
        AlgorithmLauncher.cpp
        OBJECT_LINK_LIBRARIES
        deps flags
        INTERFACE_LINK_LIBRARIES
        dmrg-algorithms-fdmrg
        dmrg-algorithms-flbit
        dmrg-algorithms-itebd
        dmrg-algorithms-idmrg
        dmrg-algorithms-xdmrg
        )

add_library(dmrg-algorithms INTERFACE)
target_link_libraries(dmrg-algorithms INTERFACE
        dmrg-algorithms-base
        dmrg-algorithms-finite
        dmrg-algorithms-infinite
        dmrg-algorithms-launcher
        dmrg-algorithms-fdmrg
        dmrg-algorithms-flbit
        dmrg-algorithms-itebd
        dmrg-algorithms-idmrg
        dmrg-algorithms-xdmrg
        )
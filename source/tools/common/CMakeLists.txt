add_subdirectory(h5)
add_subdirectory(contraction)

add_source_target(dmrg-tools-common-log TARGET_SOURCES log.cpp OBJECT_LINK_LIBRARIES dmrg-deps dmrg-flags)

add_source_target(dmrg-tools-common-prof
        TARGET_SOURCES
        prof.cpp
        OBJECT_LINK_LIBRARIES
        dmrg-deps dmrg-flags
        INTERFACE_LINK_LIBRARIES
        dmrg-tid
        dmrg-config
        dmrg-tools-common-log
        )
add_source_target(dmrg-tools-common-split
        TARGET_SOURCES
        split.cpp
        OBJECT_LINK_LIBRARIES
        dmrg-deps dmrg-flags
        INTERFACE_LINK_LIBRARIES
        dmrg-tid
        dmrg-config
        dmrg-math-svd
        dmrg-tensors-mps
        dmrg-tools-common-log
        dmrg-tools-common-contraction
        )
add_source_target(dmrg-tools-common-views
        TARGET_SOURCES
        views.cpp
        OBJECT_LINK_LIBRARIES
        dmrg-deps dmrg-flags
        INTERFACE_LINK_LIBRARIES
        dmrg-tid
        dmrg-config
        dmrg-math-eig
        dmrg-math-tenx
        dmrg-tensors-state
        dmrg-tools-common-log
        )
add_source_target(dmrg-tools-common-plot
        TARGET_SOURCES
        plot/plot.cpp
        OBJECT_LINK_LIBRARIES
        dmrg-deps dmrg-flags
        INTERFACE_LINK_LIBRARIES
        )


add_library(dmrg-tools-common INTERFACE)
target_link_libraries(dmrg-tools-common INTERFACE dmrg-deps dmrg-flags)
target_link_libraries(dmrg-tools-common INTERFACE
        dmrg-tools-common-log
        dmrg-tools-common-prof
        dmrg-tools-common-split
        dmrg-tools-common-views
        dmrg-tools-common-plot
        dmrg-tools-common-h5
        dmrg-tools-common-contraction
        )

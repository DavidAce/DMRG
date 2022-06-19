add_source_target(dmrg-tools-infinite-print
        TARGET_SOURCES
        print.cpp
        OBJECT_LINK_LIBRARIES
        deps flags
        INTERFACE_LINK_LIBRARIES
        dmrg-tensors-infinite
        )


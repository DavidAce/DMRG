add_source_target(dmrg-config
        TARGET_SOURCES
        loader.cpp
        settings.cpp
        parse.cpp
        OBJECT_LINK_LIBRARIES
        deps flags
        INTERFACE_LINK_LIBRARIES
        dmrg-tools-common-log
        )


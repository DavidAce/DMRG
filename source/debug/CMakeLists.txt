add_source_target(dmrg-debug
        TARGET_SOURCES
        stacktrace.cpp
        info.cpp
        OBJECT_LINK_LIBRARIES
        deps flags
        INTERFACE_LINK_LIBRARIES
        dmrg-tools-common-log
        )

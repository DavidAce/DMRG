add_source_target(dmrg-debug
        TARGET_SOURCES
        exceptions.cxx
        stacktrace.cpp
        info.cpp
        OBJECT_LINK_LIBRARIES
        dmrg-deps dmrg-flags
        )

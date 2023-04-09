add_source_target(dmrg-debug
        TARGET_SOURCES
        exceptions.cxx
        info.cpp
        OBJECT_LINK_LIBRARIES
        dmrg-deps dmrg-flags
        )

find_library(DW_LIBRARY NAMES dw)
if(DW_LIBRARY)
    set(BACKWARD_DW_LIBRARY ${DW_LIBRARY})
    set(BACKWARD_REDEFINE_DW BACKWARD_REDEFINE_DW=1)
endif()

add_source_target(dmrg-stacktrace
        TARGET_SOURCES
        stacktrace.cpp
        OBJECT_LINK_LIBRARIES
        dmrg-deps dmrg-flags
        PRIVATE_LINK_LIBRARIES
        Backward::Backward ${BACKWARD_DW_LIBRARY}
        COMPILE_DEFINITIONS
        ${BACKWARD_REDEFINE_DW}
        )

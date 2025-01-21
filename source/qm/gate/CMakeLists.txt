add_source_target(xdmrg++-qm-gate
        TARGET_SOURCES
        gate.cpp
        swap.cpp
        circuit.cpp
        OBJECT_LINK_LIBRARIES
        xdmrg++-deps xdmrg++-flags
        INTERFACE_LINK_LIBRARIES
        xdmrg++-math-tenx
        xdmrg++-tools-common-log
#        COMPILE_OPTIONS
#        -mlong-double-128
        )
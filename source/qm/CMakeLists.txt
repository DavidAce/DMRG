add_subdirectory(gate)
add_subdirectory(lbit)
add_subdirectory(time)
add_subdirectory(mpo)
add_subdirectory(spin)


add_library(xdmrg++-qm INTERFACE)
target_link_libraries(xdmrg++-qm INTERFACE
        xdmrg++-qm-gate
        xdmrg++-qm-spin
        xdmrg++-qm-mpo
        xdmrg++-qm-time
        xdmrg++-qm-lbit
        )

target_link_libraries(xdmrg++-qm INTERFACE xdmrg++-deps xdmrg++-flags)

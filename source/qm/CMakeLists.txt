
add_subdirectory(gate)
add_subdirectory(lbit)
add_subdirectory(time)
add_subdirectory(mpo)
add_subdirectory(spin)


add_library(dmrg-qm INTERFACE)
target_link_libraries(dmrg-qm INTERFACE
        dmrg-qm-gate
        dmrg-qm-spin
        dmrg-qm-mpo
        dmrg-qm-time
        dmrg-qm-lbit
        )

target_link_libraries(dmrg-qm INTERFACE dmrg-deps dmrg-flags)

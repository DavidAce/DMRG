
add_subdirectory(gate)
add_subdirectory(lbit)
add_subdirectory(time)
add_subdirectory(mpo)
add_subdirectory(spin)


add_library(dmrg-qm INTERFACE)
target_link_libraries(dmrg-qm INTERFACE
        $<TARGET_OBJECTS:dmrg-qm-gate>
        $<TARGET_OBJECTS:dmrg-qm-time>
        $<TARGET_OBJECTS:dmrg-qm-lbit>
        $<TARGET_OBJECTS:dmrg-qm-mpo>
        $<TARGET_OBJECTS:dmrg-qm-spin>
        )

target_link_libraries(dmrg-qm INTERFACE deps flags)

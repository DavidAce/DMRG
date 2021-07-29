add_subdirectory(common)
add_subdirectory(finite)
add_subdirectory(infinite)

add_library(dmrg-tools INTERFACE)
target_link_libraries(dmrg-tools INTERFACE dmrg-deps dmrg-flags)
target_link_libraries(dmrg-tools INTERFACE
        dmrg-tools-common
        dmrg-tools-finite
        dmrg-tools-infinite
        )


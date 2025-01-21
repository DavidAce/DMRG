add_subdirectory(common)
add_subdirectory(finite)
add_subdirectory(infinite)

add_library(xdmrg++-tools INTERFACE)
target_link_libraries(xdmrg++-tools INTERFACE xdmrg++-deps xdmrg++-flags)
target_link_libraries(xdmrg++-tools INTERFACE
        xdmrg++-tools-common
        xdmrg++-tools-finite
        xdmrg++-tools-infinite
        )


add_subdirectory(h5)
add_subdirectory(measure)
add_subdirectory(mps)
add_subdirectory(env)
add_subdirectory(opt)
add_subdirectory(print)

add_library(xdmrg++-tools-infinite INTERFACE)
target_link_libraries(xdmrg++-tools-infinite INTERFACE
        xdmrg++-tools-infinite-h5
        xdmrg++-tools-infinite-measure
        xdmrg++-tools-infinite-env
        xdmrg++-tools-infinite-opt
        xdmrg++-tools-infinite-mps
        xdmrg++-tools-infinite-print
        )
target_link_libraries(xdmrg++-tools-infinite INTERFACE xdmrg++-deps xdmrg++-flags)

add_subdirectory(h5)
add_subdirectory(measure)
add_subdirectory(mps)
add_subdirectory(env)
add_subdirectory(opt)
add_subdirectory(print)

add_library(dmrg-tools-infinite INTERFACE)
target_link_libraries(dmrg-tools-infinite INTERFACE
        dmrg-tools-infinite-h5
        dmrg-tools-infinite-measure
        dmrg-tools-infinite-env
        dmrg-tools-infinite-opt
        dmrg-tools-infinite-mps
        dmrg-tools-infinite-print
        )
target_link_libraries(dmrg-tools-infinite INTERFACE dmrg-deps dmrg-flags)

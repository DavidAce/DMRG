add_subdirectory(h5)
add_subdirectory(measure)
add_subdirectory(mps)
add_subdirectory(env)
add_subdirectory(opt)
add_subdirectory(print)

add_library(dmrg-tools-infinite INTERFACE)
target_link_libraries(dmrg-tools-infinite INTERFACE
        $<TARGET_OBJECTS:dmrg-tools-infinite-h5>
        $<TARGET_OBJECTS:dmrg-tools-infinite-measure>
        $<TARGET_OBJECTS:dmrg-tools-infinite-env>
        $<TARGET_OBJECTS:dmrg-tools-infinite-opt>
        $<TARGET_OBJECTS:dmrg-tools-infinite-mps>
        $<TARGET_OBJECTS:dmrg-tools-infinite-print>
        )
target_link_libraries(dmrg-tools-infinite INTERFACE deps flags)

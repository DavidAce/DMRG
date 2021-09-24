add_subdirectory(h5)
add_subdirectory(contraction)

add_library(dmrg-tools-common-log OBJECT)
target_sources(dmrg-tools-common-log PRIVATE log.cpp)
target_link_libraries(dmrg-tools-common-log PUBLIC deps flags)

add_library(dmrg-tools-common-prof OBJECT)
target_sources(dmrg-tools-common-prof PRIVATE prof.cpp)
target_link_libraries(dmrg-tools-common-prof PUBLIC deps flags)

add_library(dmrg-tools-common-split OBJECT)
target_sources(dmrg-tools-common-split PRIVATE split.cpp)
target_link_libraries(dmrg-tools-common-split PUBLIC deps flags)

add_library(dmrg-tools-common-views OBJECT)
target_sources(dmrg-tools-common-views PRIVATE views.cpp)
target_link_libraries(dmrg-tools-common-views PUBLIC deps flags)

target_link_precompiled_headers(dmrg-tools-common-log)
target_link_precompiled_headers(dmrg-tools-common-prof)
target_link_precompiled_headers(dmrg-tools-common-split)
target_link_precompiled_headers(dmrg-tools-common-views)


add_library(dmrg-tools-common INTERFACE)
target_link_libraries(dmrg-tools-common INTERFACE deps flags)
target_link_libraries(dmrg-tools-common INTERFACE
        $<TARGET_OBJECTS:dmrg-tools-common-log>
        $<TARGET_OBJECTS:dmrg-tools-common-prof>
        $<TARGET_OBJECTS:dmrg-tools-common-split>
        $<TARGET_OBJECTS:dmrg-tools-common-views>
        $<TARGET_OBJECTS:dmrg-tools-common-h5>
        $<TARGET_OBJECTS:dmrg-tools-common-contraction>
        )

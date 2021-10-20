add_library(dmrg-debug-o OBJECT)
target_sources(dmrg-debug-o PRIVATE stacktrace.cpp info.cpp)
target_link_libraries(dmrg-debug-o PUBLIC deps flags)
target_link_precompiled_headers(dmrg-debug-o)

add_library(dmrg-debug INTERFACE)
target_link_libraries(dmrg-debug INTERFACE $<TARGET_OBJECTS:dmrg-debug-o>)
target_link_libraries(dmrg-debug INTERFACE deps flags)

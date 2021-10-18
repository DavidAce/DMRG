add_library(dmrg-config-o OBJECT)
target_sources(dmrg-config-o PUBLIC loader.cpp settings.cpp parse.cpp)
target_link_libraries(dmrg-config-o PUBLIC deps flags)
target_link_precompiled_headers(dmrg-config-o)


add_library(dmrg-config INTERFACE)
target_link_libraries(dmrg-config INTERFACE $<TARGET_OBJECTS:dmrg-config-o>)
target_link_libraries(dmrg-config INTERFACE deps flags)

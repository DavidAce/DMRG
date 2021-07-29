add_library(dmrg-config-o OBJECT)
target_sources(dmrg-config-o PUBLIC loader.cpp settings.cpp)
target_link_libraries(dmrg-config-o PUBLIC dmrg-deps dmrg-flags)
if (DMRG_ENABLE_PCH)
    target_precompile_headers(dmrg-config-o REUSE_FROM dmrg-pch)
endif ()


add_library(dmrg-config INTERFACE)
target_link_libraries(dmrg-config INTERFACE $<TARGET_OBJECTS:dmrg-config-o>)
target_link_libraries(dmrg-config INTERFACE dmrg-deps dmrg-flags)

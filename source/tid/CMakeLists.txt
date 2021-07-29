add_library(dmrg-tid-o OBJECT)
target_sources(dmrg-tid-o PRIVATE tid.cpp ur.cpp token.cpp)
target_link_libraries(dmrg-tid-o PUBLIC dmrg-deps dmrg-flags)
if (DMRG_ENABLE_PCH)
    target_precompile_headers(dmrg-tid-o REUSE_FROM dmrg-pch)
endif ()

add_library(dmrg-tid INTERFACE)
target_link_libraries(dmrg-tid INTERFACE $<TARGET_OBJECTS:dmrg-tid-o>)
target_link_libraries(dmrg-tid INTERFACE dmrg-deps dmrg-flags)

# pronto
# ur
# yalla
# haste
# stund
# ur
# blink

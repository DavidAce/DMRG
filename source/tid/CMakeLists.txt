add_library(dmrg-tid-o OBJECT)
target_sources(dmrg-tid-o PRIVATE tid.cpp ur.cpp token.cpp)
target_link_libraries(dmrg-tid-o PUBLIC deps flags)
target_link_precompiled_headers(dmrg-tid-o)

add_library(dmrg-tid INTERFACE)
target_link_libraries(dmrg-tid INTERFACE $<TARGET_OBJECTS:dmrg-tid-o>)
target_link_libraries(dmrg-tid INTERFACE deps flags)

# pronto
# ur
# yalla
# haste
# stund
# ur
# blink

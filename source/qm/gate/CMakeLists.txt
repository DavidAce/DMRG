add_library(dmrg-qm-gate OBJECT)
target_sources(dmrg-qm-gate PRIVATE gate.cpp swap.cpp)
target_link_libraries(dmrg-qm-gate PUBLIC deps flags)
target_link_precompiled_headers(dmrg-qm-gate)

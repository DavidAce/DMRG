add_library(dmrg-qm-spin OBJECT)
target_sources(dmrg-qm-spin PRIVATE spin.cpp half.cpp one.cpp)
target_link_libraries(dmrg-qm-spin PUBLIC deps flags)
target_link_precompiled_headers(dmrg-qm-spin)

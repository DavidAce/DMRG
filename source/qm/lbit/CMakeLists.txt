add_library(dmrg-qm-lbit OBJECT)
target_sources(dmrg-qm-lbit PRIVATE lbit.cpp)
target_link_libraries(dmrg-qm-lbit PUBLIC deps flags)
target_link_precompiled_headers(dmrg-qm-lbit)

add_library(dmrg-qm-time OBJECT)
target_sources(dmrg-qm-time PRIVATE time.cpp)
target_link_libraries(dmrg-qm-time PUBLIC deps flags)
target_link_precompiled_headers(dmrg-qm-time)

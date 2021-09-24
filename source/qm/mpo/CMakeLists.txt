add_library(dmrg-qm-mpo OBJECT)
target_sources(dmrg-qm-mpo PRIVATE mpo.cpp)
target_link_libraries(dmrg-qm-mpo PUBLIC deps flags)
target_link_precompiled_headers(dmrg-qm-mpo)

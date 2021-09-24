add_library(dmrg-tools-infinite-mps OBJECT)
target_sources(dmrg-tools-infinite-mps PRIVATE mps.cpp)
target_link_libraries(dmrg-tools-infinite-mps PUBLIC deps flags)
target_link_precompiled_headers(dmrg-tools-infinite-mps)
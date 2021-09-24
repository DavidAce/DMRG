add_library(dmrg-tools-infinite-h5 OBJECT)
target_sources(dmrg-tools-infinite-h5 PUBLIC save.cpp)
target_link_libraries(dmrg-tools-infinite-h5 PUBLIC deps flags)
target_link_precompiled_headers(dmrg-tools-infinite-h5)

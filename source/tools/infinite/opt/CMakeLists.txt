add_library(dmrg-tools-infinite-opt OBJECT)
target_sources(dmrg-tools-infinite-opt PRIVATE opt.cpp)
target_link_libraries(dmrg-tools-infinite-opt PUBLIC deps flags)
target_link_precompiled_headers(dmrg-tools-infinite-opt)

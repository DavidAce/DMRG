add_library(dmrg-tools-infinite-print OBJECT)
target_sources(dmrg-tools-infinite-print PRIVATE print.cpp)
target_link_libraries(dmrg-tools-infinite-print PUBLIC deps flags)
target_link_precompiled_headers(dmrg-tools-infinite-print)

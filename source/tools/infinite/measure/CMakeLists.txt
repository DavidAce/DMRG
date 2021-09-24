add_library(dmrg-tools-infinite-measure OBJECT)
target_sources(dmrg-tools-infinite-measure PRIVATE measure.cpp measure_ham.cpp measure_mom.cpp)
target_link_libraries(dmrg-tools-infinite-measure PUBLIC deps flags)
target_link_precompiled_headers(dmrg-tools-infinite-measure)

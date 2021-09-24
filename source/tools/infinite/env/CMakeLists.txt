add_library(dmrg-tools-infinite-env OBJECT)
target_sources(dmrg-tools-infinite-env PRIVATE env.cpp)
target_link_libraries(dmrg-tools-infinite-env PUBLIC deps flags)
target_link_precompiled_headers(dmrg-tools-infinite-env)

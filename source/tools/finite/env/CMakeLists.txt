add_library(dmrg-tools-finite-env OBJECT)
target_sources(dmrg-tools-finite-env PRIVATE env.cpp)
target_link_libraries(dmrg-tools-finite-env PUBLIC deps flags)
target_link_precompiled_headers(dmrg-tools-finite-env)
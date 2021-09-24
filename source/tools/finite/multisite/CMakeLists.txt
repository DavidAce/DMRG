add_library(dmrg-tools-finite-multisite OBJECT)
target_sources(dmrg-tools-finite-multisite PRIVATE multisite.cpp)
target_link_libraries(dmrg-tools-finite-multisite PUBLIC deps flags)
target_link_precompiled_headers(dmrg-tools-finite-multisite)
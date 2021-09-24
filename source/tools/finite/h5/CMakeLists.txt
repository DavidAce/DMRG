add_library(dmrg-tools-finite-h5 OBJECT)
target_sources(dmrg-tools-finite-h5 PRIVATE load.cpp save.cpp)
target_link_libraries(dmrg-tools-finite-h5 PUBLIC deps flags)
target_link_precompiled_headers(dmrg-tools-finite-h5)
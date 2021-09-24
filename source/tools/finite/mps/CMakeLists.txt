add_library(dmrg-tools-finite-mps OBJECT)
target_sources(dmrg-tools-finite-mps PRIVATE mps.cpp mps_init_entangled_state.cpp mps_init_product_state.cpp)
target_link_libraries(dmrg-tools-finite-mps PUBLIC deps flags)
target_link_precompiled_headers(dmrg-tools-finite-mps)

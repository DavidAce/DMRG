add_library(dmrg-tools-finite-measure OBJECT)
target_sources(dmrg-tools-finite-measure PRIVATE measure.cpp measure-parity-components.cpp entropy.cpp)
target_link_libraries(dmrg-tools-finite-measure PUBLIC deps flags)
target_link_precompiled_headers(dmrg-tools-finite-measure)
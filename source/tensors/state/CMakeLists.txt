add_library(dmrg-tensors-state OBJECT)
target_sources(dmrg-tensors-state PRIVATE StateFinite.cpp StateInfinite.cpp)
target_link_libraries(dmrg-tensors-state PUBLIC deps flags)
target_link_precompiled_headers(dmrg-tensors-state)

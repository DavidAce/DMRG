add_library(dmrg-tensors-mps OBJECT)
target_sources(dmrg-tensors-mps PRIVATE MpsSite.cpp)
target_link_libraries(dmrg-tensors-mps PUBLIC deps flags)
target_link_precompiled_headers(dmrg-tensors-mps)

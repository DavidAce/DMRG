add_library(dmrg-tensors-model OBJECT)
target_sources(dmrg-tensors-model PRIVATE ModelFinite.cpp ModelInfinite.cpp)
target_link_libraries(dmrg-tensors-model PUBLIC deps flags)
target_link_precompiled_headers(dmrg-tensors-model)
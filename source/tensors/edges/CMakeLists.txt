add_library(dmrg-tensors-edges OBJECT)
target_sources(dmrg-tensors-edges PRIVATE EdgesFinite.cpp EdgesInfinite.cpp)
target_link_libraries(dmrg-tensors-edges PUBLIC deps flags)
target_link_precompiled_headers(dmrg-tensors-edges)

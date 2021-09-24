add_library(dmrg-tensors-mpo OBJECT)
target_sources(dmrg-tensors-mpo PRIVATE MpoSite.cpp MpoFactory.cpp IsingSdual.cpp IsingTfRf.cpp LBit.cpp)
target_link_libraries(dmrg-tensors-mpo PUBLIC deps flags)
target_link_precompiled_headers(dmrg-tensors-mpo)

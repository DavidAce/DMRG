add_library(dmrg-tensors-env OBJECT)
target_sources(dmrg-tensors-env PRIVATE EnvPair.cpp EnvBase.cpp EnvEne.cpp EnvVar.cpp)
target_link_libraries(dmrg-tensors-env PUBLIC deps flags)
target_link_precompiled_headers(dmrg-tensors-env)

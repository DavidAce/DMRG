add_library(dmrg-math-hash OBJECT)
target_sources(dmrg-math-hash PRIVATE hash.cpp)
target_link_libraries(dmrg-math-hash PUBLIC deps flags)
target_link_precompiled_headers(dmrg-math-hash)


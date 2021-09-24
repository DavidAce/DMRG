add_library(dmrg-math-rnd OBJECT)
target_sources(dmrg-math-rnd PRIVATE rnd.cpp)
target_link_libraries(dmrg-math-rnd PUBLIC deps flags)
target_link_precompiled_headers(dmrg-math-rnd)

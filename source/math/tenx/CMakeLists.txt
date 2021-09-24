add_library(dmrg-math-tenx OBJECT)
target_sources(dmrg-math-tenx PRIVATE omp.cpp)
target_link_libraries(dmrg-math-tenx PUBLIC deps flags)
target_link_precompiled_headers(dmrg-math-tenx)

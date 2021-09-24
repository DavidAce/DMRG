add_library(dmrg-math-linalg OBJECT)
target_sources(dmrg-math-linalg PRIVATE matrix.cpp)
target_link_libraries(dmrg-math-linalg PUBLIC deps flags)
target_link_precompiled_headers(dmrg-math-linalg)

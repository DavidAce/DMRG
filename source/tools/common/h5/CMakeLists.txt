add_library(dmrg-tools-common-h5 OBJECT)
target_sources(dmrg-tools-common-h5 PRIVATE save.cpp load.cpp resume.cpp tmp.cpp)
target_link_libraries(dmrg-tools-common-h5 PUBLIC deps flags)
target_link_precompiled_headers(dmrg-tools-common-h5)

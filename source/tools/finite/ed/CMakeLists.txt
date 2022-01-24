add_library(dmrg-tools-finite-ed OBJECT)
target_sources(dmrg-tools-finite-ed PRIVATE
        ed.cpp
        )

target_link_libraries(dmrg-tools-finite-ed PUBLIC deps flags)
target_link_precompiled_headers(dmrg-tools-finite-ed)

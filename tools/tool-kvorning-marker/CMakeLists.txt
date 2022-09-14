cmake_minimum_required(VERSION 3.15)
project(tool-kvorning-marker
        VERSION 0.0.1
        DESCRIPTION "Calculate the kvorning marker on a given MPS state"
        HOMEPAGE_URL "https://github.com/DavidAce/DMRG"
        )

add_executable(dmrg-tool-kvorning-marker main.cpp)
target_link_libraries(dmrg-tool-kvorning-marker PRIVATE dmrg-deps dmrg-flags)
target_include_directories(dmrg-tool-kvorning-marker PRIVATE ${CMAKE_SOURCE_DIR}/source)
target_link_precompiled_headers(dmrg-tool-kvorning-marker)
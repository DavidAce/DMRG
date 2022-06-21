include(../cmake/AddSourceTarget.cmake)
add_subdirectory(algorithms)
add_subdirectory(config)
add_subdirectory(debug)
add_subdirectory(io)
add_subdirectory(math)
add_subdirectory(qm)
add_subdirectory(tensors)
add_subdirectory(tid)
add_subdirectory(tools)
add_library(dmrg-main OBJECT)
target_sources(dmrg-main PRIVATE main.cpp)
target_link_libraries(dmrg-main PUBLIC dmrg-deps dmrg-flags)
target_link_libraries(dmrg-main INTERFACE
        dmrg-algorithms
        dmrg-config
        dmrg-debug
        dmrg-io
        dmrg-math
        dmrg-qm
        dmrg-tensors
        dmrg-tid
        dmrg-tools
        )

target_link_precompiled_headers(dmrg-main)
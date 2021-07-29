# This one is special, since it absolutely can't have BLAS/LAPACK headers during compile for some reason
add_library(dmrg-math-eig-arpack OBJECT)
target_sources(dmrg-math-eig-arpack PRIVATE solver_arpack.cpp)
set_source_files_properties(solver_arpack.cpp PROPERTIES SKIP_PRECOMPILE_HEADERS ON)
target_link_libraries(dmrg-math-eig-arpack PUBLIC spdlog::spdlog arpack::arpack++ dmrg-flags)

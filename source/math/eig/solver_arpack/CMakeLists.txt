add_source_target(dmrg-math-eig-arpack
        TARGET_SOURCES
        solver_arpack.cpp
        OBJECT_LINK_LIBRARIES
        deps flags
        INTERFACE_LINK_LIBRARIES
        dmrg-tid
        dmrg-math-eig-matvec
        dmrg-math-eig-o
        )
# This one is special, since it absolutely can't have BLAS/LAPACK headers during compile for some reason
set_source_files_properties(solver_arpack.cpp PROPERTIES SKIP_PRECOMPILE_HEADERS ON)




add_source_target(dmrg-math-eig-lapacke
        TARGET_SOURCES
        # Standard problems
        solver_dgeev.cpp
        solver_dsyevd.cpp
        solver_dsyevx.cpp
        solver_dsyevr.cpp
        solver_zgeev.cpp
        solver_zheevd.cpp
        solver_zheevr.cpp
        # Generalized problems
        solver_dsygvd.cpp
        solver_dsygvx.cpp

        OBJECT_LINK_LIBRARIES
        dmrg-deps dmrg-flags
        INTERFACE_LINK_LIBRARIES
        dmrg-math-eig-o
        )
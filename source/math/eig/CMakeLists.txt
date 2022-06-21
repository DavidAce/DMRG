
add_subdirectory(solver_arpack)
add_subdirectory(solver_primme)
add_subdirectory(solver_lapacke)
add_subdirectory(matvec)

add_source_target(dmrg-math-eig
        TARGET_SOURCES
        log.cpp
        settings.cpp
        solver.cpp
        solution.cpp
        OBJECT_LINK_LIBRARIES
        dmrg-deps dmrg-flags
        INTERFACE_LINK_LIBRARIES
        dmrg-tid
        dmrg-math-eig-matvec
        dmrg-math-eig-arpack
        dmrg-math-eig-lapacke
        dmrg-math-eig-primme
        )

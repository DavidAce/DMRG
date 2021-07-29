
add_library(dmrg-math-eig OBJECT)
target_sources(dmrg-math-eig PRIVATE
        log.cpp
        settings.cpp
        solver.cpp
        solution.cpp
        )

target_link_libraries(dmrg-math-eig PUBLIC dmrg-deps dmrg-flags)
if (DMRG_ENABLE_PCH)
    target_precompile_headers(dmrg-math-eig REUSE_FROM dmrg-pch)
endif ()

add_subdirectory(solver_arpack)
add_subdirectory(solver_primme)
add_subdirectory(solver_lapacke)
add_subdirectory(matvec)

target_link_libraries(dmrg-math-eig INTERFACE
        $<TARGET_OBJECTS:dmrg-math-eig-arpack>
        $<TARGET_OBJECTS:dmrg-math-eig-lapacke>
        $<TARGET_OBJECTS:dmrg-math-eig-primme>
        $<TARGET_OBJECTS:dmrg-math-eig-matvec>
        )



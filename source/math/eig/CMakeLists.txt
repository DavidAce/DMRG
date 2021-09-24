
add_library(dmrg-math-eig OBJECT)
target_sources(dmrg-math-eig PRIVATE
        log.cpp
        settings.cpp
        solver.cpp
        solution.cpp
        )

target_link_libraries(dmrg-math-eig PUBLIC deps flags)
target_link_precompiled_headers(dmrg-math-eig)


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



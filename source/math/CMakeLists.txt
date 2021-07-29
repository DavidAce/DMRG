add_subdirectory(svd)
add_subdirectory(eig)
add_subdirectory(linalg)
add_subdirectory(rnd)
add_subdirectory(hash)
add_subdirectory(tenx)
add_library(dmrg-math INTERFACE)
target_link_libraries(dmrg-math INTERFACE
        $<TARGET_OBJECTS:dmrg-math-eig>
        $<TARGET_OBJECTS:dmrg-math-eig-arpack>
        $<TARGET_OBJECTS:dmrg-math-eig-lapacke>
        $<TARGET_OBJECTS:dmrg-math-eig-primme>
        $<TARGET_OBJECTS:dmrg-math-eig-matvec>
        $<TARGET_OBJECTS:dmrg-math-svd>
        $<TARGET_OBJECTS:dmrg-math-linalg>
        $<TARGET_OBJECTS:dmrg-math-rnd>
        $<TARGET_OBJECTS:dmrg-math-hash>
        $<TARGET_OBJECTS:dmrg-math-tenx>
        )
target_link_libraries(dmrg-math INTERFACE dmrg-deps dmrg-flags)

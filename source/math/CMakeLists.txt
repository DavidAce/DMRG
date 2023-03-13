add_subdirectory(svd)
add_subdirectory(eig)
add_subdirectory(linalg)
add_subdirectory(rnd)
add_subdirectory(hash)
add_subdirectory(tenx)
add_subdirectory(fit)
add_library(dmrg-math INTERFACE)
target_link_libraries(dmrg-math INTERFACE
        dmrg-math-linalg
        dmrg-math-svd
        dmrg-math-rnd
        dmrg-math-eig
        dmrg-math-hash
        dmrg-math-tenx
        dmrg-math-fit
        )
target_link_libraries(dmrg-math INTERFACE dmrg-deps dmrg-flags)

add_source_target(dmrg-tools-common-contraction
        TARGET_SOURCES
        matrix_inverse_vector_product.cpp
        OBJECT_LINK_LIBRARIES
        dmrg-deps dmrg-flags
        INTERFACE_LINK_LIBRARIES
        dmrg-tid
        dmrg-math-tenx
        dmrg-tools-common-log
        )

if (DMRG_ENABLE_TBLIS)
    target_sources(dmrg-tools-common-contraction-o PRIVATE contraction_tblis.cpp)
    target_compile_definitions(dmrg-tools-common-contraction-o PRIVATE DMRG_ENABLE_TBLIS)
else ()
    target_sources(dmrg-tools-common-contraction-o PRIVATE contraction.cpp)
endif ()

if (DMRG_BENCH_CONTRACTION)
    target_compile_definitions(dmrg-tools-common-contraction-o PRIVATE DMRG_BENCH_CONTRACTION)
endif ()
if (DMRG_SAVE_CONTRACTION)
    target_compile_definitions(dmrg-tools-common-contraction-o PRIVATE DMRG_SAVE_CONTRACTION)
endif ()

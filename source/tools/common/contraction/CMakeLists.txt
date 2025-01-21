add_source_target(xdmrg++-tools-common-contraction
        TARGET_SOURCES
        matrix_inverse_vector_product.cpp
        OBJECT_LINK_LIBRARIES
        xdmrg++-deps xdmrg++-flags
        INTERFACE_LINK_LIBRARIES
        xdmrg++-tid
        xdmrg++-math-tenx
        xdmrg++-tools-common-log
        )

target_sources(xdmrg++-tools-common-contraction-o PRIVATE contraction_tblis.cpp)
if (DMRG_ENABLE_TBLIS)
    target_compile_definitions(xdmrg++-tools-common-contraction-o PRIVATE DMRG_ENABLE_TBLIS)
endif ()

if (DMRG_BENCH_CONTRACTION)
    target_compile_definitions(xdmrg++-tools-common-contraction-o PRIVATE DMRG_BENCH_CONTRACTION)
endif ()
if (DMRG_SAVE_CONTRACTION)
    target_compile_definitions(xdmrg++-tools-common-contraction-o PRIVATE DMRG_SAVE_CONTRACTION)
endif ()

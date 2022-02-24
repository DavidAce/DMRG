add_library(dmrg-tools-common-contraction OBJECT)
#target_sources(dmrg-tools-common-contraction PRIVATE expectation_value.cpp matrix_inverse_vector_product.cpp matrix_vector_product.cpp)
target_sources(dmrg-tools-common-contraction PRIVATE contraction.cpp matrix_inverse_vector_product_lbfgs.cpp)
target_link_libraries(dmrg-tools-common-contraction PUBLIC deps flags)
target_link_precompiled_headers(dmrg-tools-common-contraction)


if (DMRG_BENCH_CONTRACTION)
    target_compile_definitions(dmrg-tools-common-contraction PRIVATE DMRG_BENCH_CONTRACTION)
endif ()
if (DMRG_SAVE_CONTRACTION)
    target_compile_definitions(dmrg-tools-common-contraction PRIVATE DMRG_SAVE_CONTRACTION)
endif ()
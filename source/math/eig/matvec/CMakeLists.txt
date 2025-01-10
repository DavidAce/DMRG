add_source_target(dmrg-math-eig-matvec
        TARGET_SOURCES
        matvec_dense.cpp
        matvec_sparse.cpp
        matvec_mpo.cpp
        matvec_mpos.cpp
        matvec_zero.cpp
        OBJECT_LINK_LIBRARIES
        dmrg-deps dmrg-flags
        INTERFACE_LINK_LIBRARIES
        dmrg-tid
        dmrg-tensors-mps
        dmrg-tensors-mpo
        dmrg-tensors-env
        dmrg-math-svd
        dmrg-math-tenx
        dmrg-tools-common-contraction
        dmrg-math-eig-o
        )


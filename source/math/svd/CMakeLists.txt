add_source_target(dmrg-math-svd
        TARGET_SOURCES
        svd.cpp
        svd_eigen.cpp
        svd_lapacke.cpp
        svd_rsvd.cpp
        svd_save.cpp
        OBJECT_LINK_LIBRARIES
        dmrg-deps dmrg-flags
        INTERFACE_LINK_LIBRARIES
        dmrg-tid
        )


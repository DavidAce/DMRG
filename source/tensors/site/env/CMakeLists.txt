add_source_target(dmrg-tensors-env
        TARGET_SOURCES
        EnvPair.cpp
        EnvBase.cpp
        EnvEne.cpp
        EnvVar.cpp
        OBJECT_LINK_LIBRARIES
        deps flags
        INTERFACE_LINK_LIBRARIES
        dmrg-math-hash
        dmrg-tensors-mps
        dmrg-tensors-mpo
        dmrg-tensors-env
        dmrg-tools-common-log
        )
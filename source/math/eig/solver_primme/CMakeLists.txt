add_source_target(dmrg-math-eig-primme
        TARGET_SOURCES
        solver_primme.cpp
        OBJECT_LINK_LIBRARIES
        deps flags
        INTERFACE_LINK_LIBRARIES
        dmrg-tid
        dmrg-math-eig-matvec
        dmrg-math-eig-o
        )
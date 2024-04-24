add_source_target(dmrg-tensors-mpo-h5tb
                  TARGET_SOURCES
                  h5tb.cpp
                  OBJECT_LINK_LIBRARIES
                  dmrg-deps dmrg-flags
                  INTERFACE_LINK_LIBRARIES
                  dmrg-debug
                  dmrg-tools-common-log
                  )
target_include_directories(dmrg-tensors-mpo-h5tb SYSTEM INTERFACE ${CMAKE_CURRENT_LIST_DIR}) # Needed by meld
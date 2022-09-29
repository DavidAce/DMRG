add_source_target(dmrg-tools-common-h5
                  TARGET_SOURCES
                  save.cpp
                  load.cpp
                  resume.cpp
                  tmp.cpp
                  storage_info.cpp
                  OBJECT_LINK_LIBRARIES
                  dmrg-deps dmrg-flags
                  INTERFACE_LINK_LIBRARIES
                  dmrg-io
                  dmrg-tid
                  dmrg-config
                  dmrg-algorithms-status
                  dmrg-tools-common-log
                  )


add_source_target(xdmrg++-tools-common-h5
                  TARGET_SOURCES
                  save.cpp
                  load.cpp
                  resume.cpp
                  tmp.cpp
                  storage_info.cpp
                  OBJECT_LINK_LIBRARIES
                  xdmrg++-deps xdmrg++-flags
                  INTERFACE_LINK_LIBRARIES
                  xdmrg++-io
                  xdmrg++-tid
                  xdmrg++-config
                  xdmrg++-algorithms-status
                  xdmrg++-tools-common-log
                  )


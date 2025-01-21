add_source_target(xdmrg++-config
                  TARGET_SOURCES
                  loader.cpp
                  settings.cpp
                  parse.cpp
                  threading.cpp
                  OBJECT_LINK_LIBRARIES
                  xdmrg++-deps xdmrg++-flags
                  INTERFACE_LINK_LIBRARIES
                  xdmrg++-tools-common-log
                  )


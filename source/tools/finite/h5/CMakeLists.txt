add_source_target(xdmrg++-tools-finite-h5
                  TARGET_SOURCES
                  load.cpp
                  save.cpp
                  OBJECT_LINK_LIBRARIES
                  xdmrg++-deps xdmrg++-flags
                  INTERFACE_LINK_LIBRARIES
                  xdmrg++-io
                  xdmrg++-tid
                  xdmrg++-config
                  xdmrg++-algorithms-status
                  xdmrg++-math-tenx
                  xdmrg++-tensors-finite
                  xdmrg++-tools-finite-measure
                  xdmrg++-tools-finite-ops
                  xdmrg++-tools-finite-env
                  xdmrg++-tools-finite-print
                  xdmrg++-tools-common-h5
                  xdmrg++-tools-common-log
                  )


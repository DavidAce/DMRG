add_source_target(dmrg-io
                  TARGET_SOURCES
                  fmt.cxx
                  hdf5_types.cpp
                  OBJECT_LINK_LIBRARIES
                  dmrg-deps dmrg-flags
                  )

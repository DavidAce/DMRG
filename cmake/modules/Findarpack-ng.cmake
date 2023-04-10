
if(NOT TARGET arpack AND NOT TARGET ARPACK::ARPACK)
    find_package(arpackng ${arpack-ng_FIND_VERSION} CONFIG BYPASS_PROVIDER)
endif()
if(NOT TARGET arpack AND NOT TARGET arpack-ng::arpack-ng)
    find_package(arpack-ng ${arpack-ng_FIND_VERSION} CONFIG BYPASS_PROVIDER)
endif()

if(TARGET arpack-ng::arpack-ng AND NOT TARGET ARPACK::ARPACK)
    add_library(ARPACK::ARPACK ALIAS arpack-ng::arpack-ng)
endif()

if(arpack-ng_FIND_REQUIRED AND NOT TARGET ARPACK::ARPACK)
    message(FATAL_ERROR "Could NOT find arpack-ng")
endif()

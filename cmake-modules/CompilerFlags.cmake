
message(STATUS "C compiler ${CMAKE_C_COMPILER}")
message(STATUS "FC compiler ${CMAKE_Fortran_COMPILER}")
message(STATUS "CXX compiler ${CMAKE_CXX_COMPILER}")

############################################################
### Set  the same microarchitecture for c++ and OpenBLAS ###
############################################################

if(NOT DMRG_MICROARCH)
    set(DMRG_MICROARCH "native")
endif()
if(DMRG_MICROARCH)
    if (${DMRG_MICROARCH} STREQUAL "zen")
        string(TOUPPER ${DMRG_MICROARCH} OPENBLAS_MARCH)
        set(CXX_MARCH znver1)
    elseif (${DMRG_MICROARCH} STREQUAL "native")
        set(OPENBLAS_MARCH HASWELL)
        set(CXX_MARCH native)
    else()
        string(TOUPPER ${DMRG_MICROARCH} OPENBLAS_MARCH)
        string(TOLOWER ${DMRG_MICROARCH} CXX_MARCH)
    endif()
endif()


###########################################
###  Apply RELEASE/DEBUG compile flags  ###
###########################################
if(CXX_MARCH)
    message(STATUS "Using microarchitechture: ${CXX_MARCH}")
    list(APPEND CMAKE_CXX_FLAGS            -march=${CXX_MARCH} -mtune=${CXX_MARCH})
endif()
list(APPEND CMAKE_CXX_FLAGS                )
list(APPEND CMAKE_CXX_FLAGS_RELEASE        )
list(APPEND CMAKE_CXX_FLAGS_DEBUG          -fstack-protector -D_FORTIFY_SOURCE=2 -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC)
list(APPEND CMAKE_CXX_FLAGS_RELWITHDEBINFO )
list(APPEND CMAKE_CXX_FLAGS_MINSIZEREL)


string (REPLACE " " ";" CMAKE_CXX_FLAGS_LIST                "${CMAKE_CXX_FLAGS}")
string (REPLACE " " ";" CMAKE_CXX_FLAGS_RELEASE_LIST        "${CMAKE_CXX_FLAGS_RELEASE}")
string (REPLACE " " ";" CMAKE_CXX_FLAGS_DEBUG_LIST          "${CMAKE_CXX_FLAGS_DEBUG}")
string (REPLACE " " ";" CMAKE_CXX_FLAGS_RELWITHDEBINFO_LIST "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
string (REPLACE " " ";" CMAKE_CXX_FLAGS_MINSIZEREL_LIST     "${CMAKE_CXX_FLAGS_MINSIZEREL}")

list(REMOVE_DUPLICATES CMAKE_CXX_FLAGS_LIST)
list(REMOVE_DUPLICATES CMAKE_CXX_FLAGS_RELEASE_LIST)
list(REMOVE_DUPLICATES CMAKE_CXX_FLAGS_DEBUG_LIST)
list(REMOVE_DUPLICATES CMAKE_CXX_FLAGS_RELWITHDEBINFO_LIST)
list(REMOVE_DUPLICATES CMAKE_CXX_FLAGS_MINSIZEREL_LIST)

string (REPLACE ";" " " CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS_LIST}")
string (REPLACE ";" " " CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE_LIST}")
string (REPLACE ";" " " CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG_LIST}")
string (REPLACE ";" " " CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO_LIST}")
string (REPLACE ";" " " CMAKE_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_MINSIZEREL_LIST}")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
list(APPEND CMAKE_EXE_LINKER_FLAGS -fuse-ld=gold)
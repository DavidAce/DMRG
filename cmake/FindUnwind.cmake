# - Try to find libunwind
# Once done this will define
#
#  Unwind_FOUND - system has libunwind
#  unwind::unwind - cmake target for libunwind

function(find_unwind)
    if(NOT BUILD_SHARED_LIBS)
        set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX} ${CMAKE_SHARED_LIBRARY_SUFFIX})
    endif()
    find_path(Unwind_INCLUDE_DIR NAMES unwind.h libunwind.h DOC "unwind include directory")
    find_library(Unwind_LIBRARY NAMES unwind DOC "unwind library")

    if(CMAKE_SYSTEM_PROCESSOR MATCHES "^arm")
        set(Unwind_ARCH "arm")
    elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^aarch64")
        set(Unwind_ARCH "aarch64")
    elseif(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" OR
           CMAKE_SYSTEM_PROCESSOR STREQUAL "amd64" OR
           CMAKE_SYSTEM_PROCESSOR STREQUAL "corei7-64")
        set(Unwind_ARCH "x86_64")
    elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^i.86$")
        set(Unwind_ARCH "x86")
    elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^ppc64")
        set(Unwind_ARCH "ppc64")
    elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^ppc")
        set(Unwind_ARCH "ppc32")
    elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^mips")
        set(Unwind_ARCH "mips")
    elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^hppa")
        set(Unwind_ARCH "hppa")
    elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^ia64")
        set(Unwind_ARCH "ia64")
    endif(CMAKE_SYSTEM_PROCESSOR MATCHES "^arm")

    find_library(Unwind_PLATFORM_LIBRARY NAMES "unwind-${Unwind_ARCH}"
                 DOC "unwind library platform")

    mark_as_advanced(Unwind_INCLUDE_DIR Unwind_LIBRARY Unwind_PLATFORM_LIBRARY)

    # Extract version information
    if(Unwind_LIBRARY)
        set(_Unwind_VERSION_HEADER ${Unwind_INCLUDE_DIR}/libunwind-common.h)

        if(EXISTS ${_Unwind_VERSION_HEADER})
            file(READ ${_Unwind_VERSION_HEADER} _Unwind_VERSION_CONTENTS)

            string(REGEX REPLACE ".*#define UNW_VERSION_MAJOR[ \t]+([0-9]+).*" "\\1"
                   Unwind_VERSION_MAJOR "${_Unwind_VERSION_CONTENTS}")
            string(REGEX REPLACE ".*#define UNW_VERSION_MINOR[ \t]+([0-9]+).*" "\\1"
                   Unwind_VERSION_MINOR "${_Unwind_VERSION_CONTENTS}")
            string(REGEX REPLACE ".*#define UNW_VERSION_EXTRA[ \t]+([0-9]+).*" "\\1"
                   Unwind_VERSION_PATCH "${_Unwind_VERSION_CONTENTS}")

            set(Unwind_VERSION ${Unwind_VERSION_MAJOR}.${Unwind_VERSION_MINOR})

            if(CMAKE_MATCH_0)
                # Third version component may be empty
                set(Unwind_VERSION ${Unwind_VERSION}.${Unwind_VERSION_PATCH})
                set(Unwind_VERSION_COMPONENTS 3)
            else(CMAKE_MATCH_0)
                set(Unwind_VERSION_COMPONENTS 2)
            endif(CMAKE_MATCH_0)
        endif()
    endif()
endfunction()

# handle the QUIETLY and REQUIRED arguments and set Unwind_FOUND to TRUE
# if all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Unwind REQUIRED_VARS Unwind_INCLUDE_DIR
                                  Unwind_LIBRARY Unwind_PLATFORM_LIBRARY VERSION_VAR Unwind_VERSION)

if(Unwind_FOUND)
    if(NOT TARGET unwind::unwind)
        add_library(unwind::unwind INTERFACE IMPORTED)
        target_include_directories(unwind::unwind SYSTEM INTERFACE ${Unwind_INCLUDE_DIR})
        target_link_libraries(unwind::unwind INTERFACE ${Unwind_LIBRARY} ${Unwind_PLATFORM_LIBRARY})
        set_target_properties(unwind::unwind PROPERTIES IMPORTED_CONFIGURATIONS RELEASE)
    endif()
endif()

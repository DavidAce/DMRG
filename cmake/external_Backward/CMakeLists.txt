cmake_minimum_required(VERSION 3.15)
project(external-Backward)
if(INIT_CACHE_FILE)
    set(INIT_CACHE_ARG -C${INIT_CACHE_FILE})
endif()
include(GNUInstallDirs)
include(ExternalProject)

option(STACK_WALKING_UNWIND             "" ON)   # Better than BACKTRACE
option(STACK_WALKING_BACKTRACE          "" OFF)
option(STACK_WALKING_LIBUNWIND          "" OFF)
option(STACK_DETAILS_AUTO_DETECT        "" ON)
option(STACK_DETAILS_BACKTRACE_SYMBOL   "" OFF) # Last resort
option(STACK_DETAILS_DW                 "" ON)
option(STACK_DETAILS_BFD                "" ON)
option(STACK_DETAILS_DWARF              "" ON)
set(BUILD_SHARED_LIBS OFF CACHE INTERNAL "")

ExternalProject_Add(external_Backward
        URL https://github.com/bombela/backward-cpp/archive/refs/tags/v1.6.tar.gz
        URL_MD5 0facf6e0fb35ed0f3cd069424a1dc79a
        BUILD_ALWAYS TRUE
        PREFIX          ${CMAKE_BINARY_DIR}
        INSTALL_DIR     ${CMAKE_INSTALL_PREFIX}
        #Give Ninja access to the terminal.
        USES_TERMINAL_DOWNLOAD TRUE
        USES_TERMINAL_UPDATE TRUE
        USES_TERMINAL_CONFIGURE TRUE
        USES_TERMINAL_BUILD TRUE
        USES_TERMINAL_INSTALL TRUE
        USES_TERMINAL_TEST TRUE

        CMAKE_ARGS
        ${INIT_CACHE_ARG}
        -DBUILD_SHARED_LIBS:BOOL=OFF
        -DCMAKE_FIND_LIBRARY_SUFFIXES:STRING=${CMAKE_STATIC_LIBRARY_SUFFIX}
        -DCMAKE_POLICY_DEFAULT_CMP0074=NEW
        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
        -DCMAKE_INSTALL_RPATH:PATH=<INSTALL_DIR>/${CMAKE_INSTALL_LIBDIR}
        -DCMAKE_SKIP_BUILD_RPATH:BOOL=OFF
        -DCMAKE_BUILD_WITH_INSTALL_RPATH:BOOL=OFF
        -DCMAKE_INSTALL_RPATH_USE_LINK_PATH:BOOL=ON
        -DCMAKE_INSTALL_MESSAGE=LAZY
        # Backward flags
        -DBACKWARD_SHARED:BOOL=OFF
        -DSTACK_WALKING_UNWIND:BOOL=${STACK_WALKING_UNWIND}
        -DSTACK_WALKING_BACKTRACE:BOOL=FALSE
        -DSTACK_WALKING_LIBUNWIND:BOOL=FALSE
        -DSTACK_DETAILS_AUTO_DETECT:BOOL=${STACK_DETAILS_AUTO_DETECT}
        -DSTACK_DETAILS_BACKTRACE_SYMBOL:BOOL=FALSE
        -DSTACK_DETAILS_DW:BOOL=${STACK_DETAILS_DW}                  # I think this is newer than DWARF and more efficient https://github.com/stephenrkell/libdwarfpp/issues/12
        -DSTACK_DETAILS_BFD:BOOL=${DSTACK_DETAILS_BFD}
        -DSTACK_DETAILS_DWARF:BOOL=${STACK_DETAILS_DWARF}
        )


# Helpful comment from https://github.com/bombela/backward-cpp/issues/124#issuecomment-432271062

#bombela commented on Oct 23, 2018
#
#Sorry, my reply was really confusing. I meant that DW and DWARF are the two best options for reading debug information. Pick the one easiest to install on your system. Let me expand.
#
#There is two important dependencies that are needed for detailed stacktrace.
#The first, is the ability to walk the call stack. The second is the ability to decode debug information from the executable.
#
#I will be assuming Linux here.
#Walking the stack
#
#For walking the stack you have the choice between UNWIND and BACKTRACE (see backward.hpp#L85). UNWIND is the best possible choice. BACKTRACE doesn't report if the address points before the instruction that failed (cpu exception) or after it (signal/interruption).
#Reading debug information
#
#For decoding debug information in the DWARF format, the choice is greater (see backward.hpp#L115). The two libraries that offer the best support and most details are DW and DWARF. BFD follows. And BACKTRACE_SYMBOL comes last. If you paid attention, you will notice that BACKTRACE_SYMBOL is the same library as BACKTRACE mentioned above. It offers the bare minimum, but sometimes that's all you have.
#Conclusion
#
#UNWIND + (DW or DWARF) will give you perfectly accurate stack traces.